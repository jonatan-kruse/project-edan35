/*
 *  main.cpp
 *  swTracer
 *
 *  Created by Michael Doggett on 2021-09-23.
 *  Copyright (c) 2021 Michael Doggett
 */
#define _USE_MATH_DEFINES
#include <cfloat>
#include <cmath>
#include <ctime>
#include <iostream>
#include <random>
#include <sstream>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "swCamera.h"
#include "swIntersection.h"
#include "swMaterial.h"
#include "swRay.h"
#include "swScene.h"
#include "swSphere.h"
#include "swVec3.h"

using namespace sw;
using namespace std;

inline float clamp(float x, float min, float max) {
    if (x < min) return min;
    if (x > max) return max;
    return x;
}

float uniform() {
    // Will be used to obtain a seed for the random number engine
    static std::random_device rd;
    // Standard mersenne_twister_engine seeded with rd()
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<float> dis(0.0f, 1.0f);
    return dis(gen);
}

void writeColor(int index, Vec3 p, uint8_t *pixels) {
    // gamma correct for gamma=2.2, x^(1/gamma), more see :
    // https://www.geeks3d.com/20101001/tutorial-gamma-correction-a-story-of-linearity/
    for (int n = 0; n < 3; n++) {
        p.m[n] = pow(p.m[n], 1.0f / 2.2f);
        pixels[index + n] = (uint8_t)(256 * clamp(p.m[n], 0.0f, 0.999f));
    }
}

Color traceRay(const Ray &r, Scene scene, int depth) {
    Color c, directColor, reflectedColor, refractedColor;
    if (depth < 0) return c;

    Intersection hit, shadow;
    if (!scene.intersect(r, hit)) return Color(0.0f, 0.0f, 0.0f); // Background color

    const Vec3 lightPos(0.0f, 30.0f, -5.0f);
    Vec3 lightDir = lightPos - hit.position;
    lightDir.normalize();
    Vec3 normal = hit.normal.normalize();

    Ray shadowRay = hit.getShadowRay(lightPos);
    if (scene.intersect(shadowRay, shadow, true)) {
        directColor = Color(0.0f, 0.0f, 0.0f);
    } else {
        directColor = hit.material.color * (lightDir * normal);
    }

    Ray reflectedRay = hit.getReflectedRay();
    reflectedColor = traceRay(reflectedRay, scene, depth - 1);

    Ray refractedRay = hit.getRefractedRay();
    refractedColor = traceRay(refractedRay, scene, depth - 1);

    c = (1 - hit.material.reflectivity - hit.material.transparency) * directColor +
        hit.material.reflectivity * reflectedColor + hit.material.transparency * refractedColor;

    return c;
}

void sphere_collision(Sphere &s1, Sphere &s2) {
    Vec3 dc1 = s1.center - s2.center;
    Vec3 dv1 = s1.velocity - s2.velocity;
    Vec3 v1 = ((dv1 * dc1) / (dc1 * dc1)) * dc1;

    Vec3 dc2 = s2.center - s1.center;
    Vec3 dv2 = s2.velocity - s1.velocity;
    Vec3 v2 = ((dv2 * dc2) / (dc2 * dc2)) * dc2;

    s1.set_velocity(s1.velocity - v1);
    s2.set_velocity(s2.velocity - v2);

    // move them out from each other
    float penetration = s1.radius + s2.radius - dc1.norm();
    s1.center = s1.center + 0.5f * penetration * dc1.normalize();
    s2.center = s2.center + 0.5f * penetration * dc2.normalize();
}

int main() {
    const int imageWidth = 1080;
    const int imageHeight = imageWidth;
    const int numChannels = 3;
    int depth = 3;
    int ssGridSize = 3; // TODO set to higher when rendering for real
    float ssStep = 1.0f / ssGridSize;
    float animation_time = 10.0f;
    float fps = 34.0f;
    float physics_update_per_frame = 5.0f;
    float dt = 1.0f / (physics_update_per_frame * fps);

    // Define vertices for Cornell box
    Vec3 vertices[] = {
      Vec3(-20.0f, 0.0f, 50.0f),  Vec3(20.0f, 0.0f, 50.0f),    Vec3(20.0f, 0.0f, -50.0f),   // Floor 1
      Vec3(-20.0f, 0.0f, 50.0f),  Vec3(20.0f, 0.0f, -50.0f),   Vec3(-20.0f, 0.0f, -50.0f),  // Floor 2
      Vec3(-20.0f, 0.0f, -50.0f), Vec3(20.0f, 0.0f, -50.0f),   Vec3(20.0f, 40.0f, -50.0f),  // Back wall 1
      Vec3(-20.0f, 0.0f, -50.0f), Vec3(20.0f, 40.0f, -50.0f),  Vec3(-20.0f, 40.0f, -50.0f), // Back wall 2
      Vec3(-20.0f, 40.0f, 50.0f), Vec3(-20.0f, 40.0f, -50.0f), Vec3(20.0f, 40.0f, 50.0f),   // Ceiling 1
      Vec3(20.0f, 40.0f, 50.0f),  Vec3(-20.0f, 40.0f, -50.0f), Vec3(20.0f, 40.0f, -50.0f),  // Ceiling 2
      Vec3(-20.0f, 0.0f, 50.0f),  Vec3(-20.0f, 40.0f, -50.0f), Vec3(-20.0f, 40.0f, 50.0f),  // Red wall 1
      Vec3(-20.0f, 0.0f, 50.0f),  Vec3(-20.0f, 0.0f, -50.0f),  Vec3(-20.0f, 40.0f, -50.0f), // Red wall 2
      Vec3(20.0f, 0.0f, 50.0f),   Vec3(20.0f, 40.0f, -50.0f),  Vec3(20.0f, 40.0f, 50.0f),   // Green wall 1
      Vec3(20.0f, 0.0f, 50.0f),   Vec3(20.0f, 0.0f, -50.0f),   Vec3(20.0f, 40.0f, -50.0f)   // Green wall 2
    };

    vector<Vec3> floor = {Vec3(-20.0f, 0.0f, 50.0f), Vec3(20.0f, 0.0f, 50.0f), Vec3(20.0f, 0.0f, -50.0f),
                          Vec3(-20.0f, 0.0f, -50.0f)};
    vector<Vec3> backWall = {Vec3(-20.0f, 0.0f, -50.0f), Vec3(20.0f, 0.0f, -50.0f), Vec3(20.0f, 40.0f, -50.0f),
                             Vec3(-20.0f, 40.0f, -50.0f)};
    vector<Vec3> ceiling = {Vec3(-20.0f, 40.0f, 50.0f), Vec3(-20.0f, 40.0f, -50.0f), Vec3(20.0f, 40.0f, 50.0f),
                            Vec3(20.0f, 40.0f, -50.0f)};
    vector<Vec3> redWall = {Vec3(-20.0f, 0.0f, 50.0f), Vec3(-20.0f, 40.0f, -50.0f), Vec3(-20.0f, 40.0f, 50.0f),
                            Vec3(-20.0f, 0.0f, -50.0f)};
    vector<Vec3> greenWall = {Vec3(20.0f, 0.0f, 50.0f), Vec3(20.0f, 40.0f, 50.0f), Vec3(20.0f, 40.0f, -50.0f),
                              Vec3(20.0f, 0.0f, -50.0f)};
    vector<vector<Vec3>> walls = {floor, backWall, ceiling, redWall, greenWall};

    // Define materials
    Material whiteDiffuse = Material(Color(0.9f, 0.9f, 0.9f), 0.0f, 0.0f, 1.0f);
    Material greenDiffuse = Material(Color(0.1f, 0.6f, 0.1f), 0.0f, 0.0f, 1.0f);
    Material redDiffuse = Material(Color(1.0f, 0.1f, 0.1f), 0.0f, 0.0f, 1.0f);
    Material blueDiffuse = Material(Color(0.0f, 0.2f, 0.9f), 0.0f, 0.0f, 1.0f);
    Material yellowReflective = Material(Color(1.0f, 0.6f, 0.1f), 0.2f, 0.0f, 1.0f);
    Material transparent = Material(Color(1.0f, 1.0f, 1.0f), 0.2f, 0.8f, 1.3f);

    Sphere s1 = Sphere(Vec3(-7.0f, 3.0f, -20.0f), 3.0f, greenDiffuse);
    Sphere s2 = Sphere(Vec3(0.0f, 3.0f, -20.0f), 3.0f, blueDiffuse);
    Sphere s3 = Sphere(Vec3(7.0f, 3.0f, -20.0f), 3.0f, redDiffuse);

    Sphere reflectiveSphere1 = Sphere(Vec3(7.0f, 3.0f, 0.0f), 3.0f, yellowReflective);
    Sphere reflectiveSphere2 = Sphere(Vec3(9.0f, 10.0f, 0.0f), 3.0f, yellowReflective);

    Sphere transparentSphere1 = Sphere(Vec3(-7.0f, 3.0f, 0.0f), 3.0f, transparent);
    Sphere transparentSphere2 = Sphere(Vec3(-9.0f, 10.0f, 0.0f), 3.0f, transparent);

    s1.set_velocity(Vec3(3.0f, 30.0f, 15.0f));
    s2.set_velocity(Vec3(0.0f, 50.0f, 0.0f));
    s3.set_velocity(Vec3(15.0f, 30.0f, 0.0f));

    transparentSphere1.set_velocity(Vec3(-5.0f, 30.0f, -15.0f));
    transparentSphere2.set_velocity(Vec3(0.0f, 50.0f, -6.0f));

    reflectiveSphere1.set_velocity(Vec3(0.0f, 30.0f, -15.0f));
    reflectiveSphere2.set_velocity(Vec3(5.0f, 50.0f, -5.0f));

    vector<Sphere> spheres = {s1, s2, s3, reflectiveSphere1, reflectiveSphere2, transparentSphere1, transparentSphere2};

    float time = 0.0f;
    int physics_updates = 0;
    while (time < animation_time) {
        // Setup scene
        Scene scene;

        scene.push(Triangle(&vertices[0], whiteDiffuse)); // Floor 1
        scene.push(Triangle(&vertices[3], whiteDiffuse)); // Floor 2

        scene.push(Triangle(&vertices[6], whiteDiffuse));  // Back wall 1
        scene.push(Triangle(&vertices[9], whiteDiffuse));  // Back wall 2
        scene.push(Triangle(&vertices[12], whiteDiffuse)); // Ceiling 1
        scene.push(Triangle(&vertices[15], whiteDiffuse)); // Ceiling 2
        scene.push(Triangle(&vertices[18], redDiffuse));   // Red wall 1
        scene.push(Triangle(&vertices[21], redDiffuse));   // Red wall 2
        scene.push(Triangle(&vertices[24], greenDiffuse)); // Green wall 1
        scene.push(Triangle(&vertices[27], greenDiffuse)); // Green wall 2

        for (int i = 0; i < spheres.size(); i++) {
            Sphere &s = spheres[i];
            for (vector<Vec3> &wall : walls) {
                Vec3 normal = ((wall[1] - wall[0]) % (wall[2] - wall[0])).normalize();
                float distance_to_wall = abs((wall[0] - s.center) * normal);
                if (distance_to_wall < s.radius) {
                    s.set_velocity(s.velocity - 2.0f * (s.velocity * normal) * normal);
                    float penetration = s.radius - distance_to_wall;
                    s.center = s.center + 2.0f * penetration * normal;
                }
            }
            for (int j = i + 1; j < spheres.size(); j++) {
                Sphere &s2 = spheres[j];
                float distance = (s.center - s2.center).norm();
                if (distance < s.radius + s2.radius) {
                    sphere_collision(s, s2);
                }
            }

            // gravity (always applied)
            s.set_velocity(s.velocity + Vec3(0.0f, -dt * 98.0f, 0.0f));

            scene.push(s);
        }

        // Setup camera
        Vec3 eye(0.0f, 10.0f, 30.0f);
        Vec3 lookAt(0.0f, 10.0f, -5.0f);
        Vec3 up(0.0f, 1.0f, 0.0f);
        Camera camera(eye, lookAt, up, 52.0f, (float)imageWidth / (float)imageHeight);
        camera.setup(imageWidth, imageHeight);

        // Ray trace pixels
        if (physics_updates < physics_update_per_frame) {
            physics_updates++;
        } else {
            uint8_t *pixels = new uint8_t[imageWidth * imageHeight * numChannels];
            std::cout << "Rendering... " << time << " s\n";
            for (int j = 0; j < imageHeight; ++j) {
                for (int i = 0; i < imageWidth; ++i) {
                    Color pixel;
                    // Get center of pixel coordinate
                    float cx = ((float)i) + 0.5f;
                    float cy = ((float)j) + 0.5f;

                    for (int k = 0; k < ssGridSize; ++k) {
                        for (int l = 0; l < ssGridSize; ++l) {
                            cx = ((float)i) + ssStep * uniform() + ssStep * k;
                            cy = ((float)j) + ssStep * uniform() + ssStep * l;

                            // Get a ray and trace it
                            Ray r = camera.getRay(cx, cy);
                            pixel = pixel + traceRay(r, scene, depth) * (1.0f / (ssGridSize * ssGridSize));
                        }
                    }

                    // Write pixel value to image
                    writeColor((j * imageWidth + i) * numChannels, pixel, pixels);
                }
            }
            char *filename = new char[100];
            sprintf(filename, "output%04d.png", (int)(time * 1000));
            stbi_write_png(filename, imageWidth, imageHeight, numChannels, pixels, imageWidth * numChannels);
            delete[] pixels;
            physics_updates = 0;
        }

        for (Sphere &s : spheres) {
            s.tick(dt);
        }

        time += dt;
    }
}