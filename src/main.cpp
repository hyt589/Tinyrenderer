#include <iostream>
#include <tgaimage.h>
#include <limits>
#include "model.h"
#include "geometry.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor green = TGAColor(0, 255, 0, 255);

const int width = 800;
const int height = 800;

float *myzbuffer = new float[width * height];

void line(Vec2i t0, Vec2i t1, TGAImage &image, TGAColor color)
{

    bool steep = false;
    if (std::abs(t0.x - t1.x) < std::abs(t0.y - t1.y))
    {
        steep = true;
        std::swap(t0.x, t0.y);
        std::swap(t1.x, t1.y);
    }
    if (t0.x > t1.x)
    {
        std::swap(t0.x, t1.x);
        std::swap(t0.y, t1.y);
    }

    int y = t0.y;
    int dx = std::abs(t0.x - t1.x);
    int dy = std::abs(t0.y - t1.y);
    int derror = std::abs(t0.y - t1.y) * 2;
    int error = 0;
    for (int x = t0.x; x <= t1.x; x++)
    {
        if (steep)
        {
            image.set(y, x, color);
        }
        else
        {
            image.set(x, y, color);
        }
        error += derror;
        if (error > dx)
        {
            y += (y1 > y0 ? 1 : -1); //add or minus 1 depending on slope
            error -= 2 * dx;
        }
    }
}

void oldSchoolTriangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color)
{

    if (t0.y > t1.y)
        std::swap(t0, t1);
    if (t0.y > t2.y)
        std::swap(t0, t2);
    if (t1.y > t2.y)
        std::swap(t1, t2);
    int total_height = t2.y - t0.y;
    for (int i = 0; i < total_height; i++)
    {
        bool second_half = i > t1.y - t0.y || t1.y == t0.y;
        int segment_height = second_half ? t2.y - t1.y : t1.y - t0.y;
        float alpha = (float)i / total_height;
        float beta = (float)(i - (second_half ? t1.y - t0.y : 0)) / segment_height; // be careful: with above conditions no division by zero here
        Vec2i A = t0 + (t2 - t0) * alpha;
        Vec2i B = second_half ? t1 + (t2 - t1) * beta : t0 + (t1 - t0) * beta;
        if (A.x > B.x)
            std::swap(A, B);
        for (int j = A.x; j <= B.x; j++)
        {
            image.set(j, t0.y + i, color); // attention, due to int casts t0.y+i != A.y
        }
    }
}

template <typename T>
Vec3<T> cross(Vec3<T> a, Vec3<T> b)
{
    Vec3<T> c(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
    return c;
}

Vec3f barycentric(Vec2i pts[3], Vec2i p)
{
    Vec3f u = cross(Vec3f(pts[2].x - pts[0].x, pts[1].x - pts[0].x, pts[0].x - p.x),
                    Vec3f(pts[2].y - pts[0].y, pts[1].y - pts[0].y, pts[0].y - p.y));
    /* `pts` and `P` has integer value as coordinates
       so `abs(u[2])` < 1 means `u[2]` is 0, that means
       triangle is degenerate, in this case return something with negative coordinates */

    if (std::abs(u.z) < 1)
        return Vec3f(-1, 1, 1);
    return Vec3f(1.f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z);
}

Vec3f barycentric(Vec3f A, Vec3f B, Vec3f C, Vec3f P)
{
    Vec3f s[2];
    for (int i = 2; i--;)
    {
        s[i].raw[0] = C.raw[i] - A.raw[i];
        s[i].raw[1] = B.raw[i] - A.raw[i];
        s[i].raw[2] = A.raw[i] - P.raw[i];
    }
    Vec3f u = cross(s[0], s[1]);
    if (std::abs(u.raw[2]) > 1e-2) // dont forget that u[2] is integer. If it is zero then triangle ABC is degenerate
        return Vec3f(1.f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z);
    return Vec3f(-1, 1, 1); // in this case generate negative coordinates, it will be thrown away by the rasterizator
}

/*
For each pixel we compute its barycentric coordinates. If it has at least one negative component, then the pixel is outside of the triangle.
*/

void triangle(Vec2i pts[3], TGAImage &image, TGAColor color)
{
    Vec2i bboxmin(image.get_width() - 1, image.get_height() - 1);
    Vec2i bboxmax(0, 0);
    Vec2i clamp(image.get_width() - 1, image.get_height() - 1);
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            bboxmin.raw[j] = std::max(0, std::min(bboxmin.raw[j], pts[i].raw[j]));
            bboxmax.raw[j] = std::min(clamp.raw[j], std::max(bboxmax.raw[j], pts[i].raw[j]));
        }
    }
    Vec2i P;
    for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++)
    {
        for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++)
        {
            Vec3f bc_screen = barycentric(pts, P);
            if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0)
                continue;
            image.set(P.x, P.y, color);
        }
    }
}

void triangle(Vec3f pts[3], float *zbuffer, TGAImage &image, TGAColor color)
{
    Vec2f bboxmin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
    Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    Vec2f clamp(image.get_width() - 1, image.get_height() - 1);
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            bboxmin.raw[j] = std::max(0.f, std::min(bboxmin.raw[j], pts[i].raw[j]));
            bboxmax.raw[j] = std::min(clamp.raw[j], std::max(bboxmax.raw[j], pts[i].raw[j]));
        }
    }
    Vec3f P;
    for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++)
    {
        for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++)
        {
            Vec3f bc_screen = barycentric(pts[0], pts[1], pts[2], P);
            if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0)
                continue;
            P.z = 0;
            for (int i = 0; i < 3; i++)
                P.z += pts[i].raw[2] * bc_screen.raw[i];
            if (zbuffer[int(P.x + P.y * width)] < P.z)
            {
                zbuffer[int(P.x + P.y * width)] = P.z;
                image.set(P.x, P.y, color);
            }
        }
    }
}


Vec2i xyFrom(Vec3f v){
    return Vec2i(int(v.x), int(v.y));
}

void triangle(Vec3f pts[3], float *zbuffer, TGAImage &image, Vec2i text_mapping[3], TGAImage &texture, float intensity)
{
    Vec2f bboxmin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
    Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    Vec2f clamp(image.get_width() - 1, image.get_height() - 1);
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            bboxmin.raw[j] = std::max(0.f, std::min(bboxmin.raw[j], pts[i].raw[j]));
            bboxmax.raw[j] = std::min(clamp.raw[j], std::max(bboxmax.raw[j], pts[i].raw[j]));
        }
    }
    Vec3f P;
    for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++)
    {
        for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++)
        {
            Vec3f bc_screen = barycentric(pts[0], pts[1], pts[2], P);
            if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0)
                continue;
            P.z = 0;
            for (int i = 0; i < 3; i++)
                P.z += pts[i].raw[2] * bc_screen.raw[i];
            if (zbuffer[int(P.x + P.y * width)] < P.z)
            {
                zbuffer[int(P.x + P.y * width)] = P.z;
                Vec2i ab_t = text_mapping[1] - text_mapping[0];
                Vec2i ac_t = text_mapping[2] - text_mapping[0];
                Vec2i color_t = ab_t * bc_screen.y + ac_t * bc_screen.z + text_mapping[0];
                // std::cout << color_t << std::endl;
                image.set(P.x, P.y, TGAColor(texture.get(color_t.x, color_t.y).r * intensity, texture.get(color_t.x, color_t.y).g * intensity, texture.get(color_t.x, color_t.y).b * intensity, 255));
            }
        }
    }
}

Vec3f world2screen(Vec3f v)
{
    return Vec3f(int((v.x + 1.) * width / 2. + .5), int((v.y + 1.) * height / 2. + .5), v.z);
}

#include <fstream>
#include <regex>
#include <tuple>

std::tuple<std::vector<Vec2i>, std::vector<std::vector<int>>> textCoords()
{
    std::ifstream file("obj/african_head.obj");
    std::vector<Vec2i> vts;
    std::vector<std::vector<int>> faces;
    if (file.is_open())
    {
        std::string line;
        while (std::getline(file, line))
        {
            // using printf() in all tests for consistency
            std::regex vt("(vt)(.*)");
            std::regex f("f(.*)");
            if (std::regex_match(line, vt))
            {
                Vec2f posf(std::stof(line.substr(5, 5)), std::stof(line.substr(11, 5)));
                Vec2i posi(int(posf.x * 1024), int(posf.y * 1024));
                vts.push_back(posi);
            } else if (std::regex_match(line, f))
            {
                std::regex mid("/(\\d+)/");
                std::vector<int> f;

                std::sregex_iterator iter(line.begin(), line.end(), mid);
                std::sregex_iterator end;
                while (iter != end)
                {
                    f.push_back(std::stoi((*iter)[1])-1);            
                    iter++;
                }
                faces.push_back(f);
            }
        }
        file.close();
    }

    std::cout << "# vt: " << vts.size() << " # f: " << faces.size() << std::endl;
    // for (size_t i = 0; i < vts.size(); i++)
    // {
    //     std::cout << vts[i] << std::endl;
    // }
    

    // for (size_t i = 0; i < faces.size(); i++)
    // {
    //     for (size_t j = 0; j < faces[i].size(); j++)
    //     {
    //         std::cout << faces[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    

    return {vts, faces};
}

int main(int argc, char **argv)
{

    auto [text_coords, text_faces] = textCoords();

    TGAImage texture;
    bool texture_loaded = texture.read_tga_file("obj/african_head_diffuse.tga");
    texture.flip_vertically();

    for (size_t i = 0; i < width * height; i++)
    {
        myzbuffer[i] = -std::numeric_limits<float>::max();
    }

    Model model("obj/african_head.obj");
    TGAImage image(width, height, TGAImage::RGB);

    for (int i = 0; i < model.nfaces(); i++)
    {
        std::vector<int> face = model.face(i);
        Vec3f screen_coords_3f[3];
        Vec3f world_coords[3];
        Vec2i text_mapping[3];
        for (int j = 0; j < 3; j++)
        {
            Vec3f v = model.vert(face[j]);
            world_coords[j] = v;
            screen_coords_3f[j] = world2screen(v);
            text_mapping[j] = text_coords[text_faces[i][j]];
            // std::cout << text_faces[i][j] << "/" << face[j] << " ";
        }
        // std::cout << std::endl;
        Vec3f n = (world_coords[2] - world_coords[0]) ^ (world_coords[1] - world_coords[0]);
        n.normalize();
        Vec3f light_dir(1.0, -1.0, -1.0);
        light_dir.normalize();
        float intensity = n * light_dir;
        if (intensity > 0) // if intensity < 0 => light from beghind polygon, discard this triangle -- "back-face culling"
        {
            triangle(screen_coords_3f, myzbuffer, image, text_mapping, texture, intensity);
        }
    }

    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output.tga");
    return 0;
}

void modernTriangle();