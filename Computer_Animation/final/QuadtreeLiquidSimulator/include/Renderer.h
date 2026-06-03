#include "Base.h"
#include "QT.h"
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <raylib.h>
#include <typeinfo>

class FluidRenderer {
  private:
    const BaseSimulator &sim;
    int screenWidth, screenHeight;
    Image image;
    Texture2D texture;
    Color *pixels;

  public:
    FluidRenderer(const BaseSimulator &sim, int screenW, int screenH)
        : sim(sim), screenWidth(screenW), screenHeight(screenH) {

        // 建立與模擬網格大小相同的 Image
        pixels = new Color[sim.getWidth() * sim.getHeight()];
        image  = {
            pixels,
            sim.getWidth(),
            sim.getHeight(),
            1,
            PIXELFORMAT_UNCOMPRESSED_R8G8B8A8
        };
        texture = LoadTextureFromImage(image);
    }

    ~FluidRenderer() {
        UnloadTexture(texture);
        delete[] pixels;
    }

    void draw(
        bool showGrid = false, bool showParticle = false, bool showPhi = false
    ) {
        int nx = sim.getWidth();
        int ny = sim.getHeight();

        // 在 CPU 端快速填充像素
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < ny; i++) {
                int idx = i + j * nx;
                if (sim.is_water(i, j)) {
                    pixels[idx] = Color{0, 50, 100, 255}; // 你的科技螢光藍
                } else {
                    pixels[idx] = BLANK; // 透明或黑色
                }
            }
        }

        // 一次性將資料送給 GPU
        UpdateTexture(texture, pixels);

        // 放大畫回螢幕上
        Rectangle source = {0, 0, (float)nx, (float)ny};
        Rectangle dest   = {0, 0, (float)screenWidth, (float)screenHeight};
        DrawTexturePro(texture, source, dest, {0, 0}, 0.0f, WHITE);

        float scaleX = (float)screenWidth / nx;
        float scaleY = (float)screenHeight / ny;

        if (showPhi) {
            try {
                const QTSimulator &qt_sim =
                    dynamic_cast<const QTSimulator &>(sim);
                for (int j = 0; j < ny; j++) {
                    for (int i = 0; i < nx; i++) {
                        int idx = i + j * nx;
                        float phi =
                            qt_sim.getNodeAt((float)i + 0.5f, (float)j + 0.5f)
                                .phi;

                        Color color;
                        if (phi <= 0.f)
                            color = Color{
                                0, (unsigned char)(20.f * std::abs(phi)), 0, 100
                            };
                        else
                            color = Color{
                                (unsigned char)(20.f * std::abs(phi)), 0, 0, 100
                            };

                        pixels[idx] = color;
                    }
                }
                UpdateTexture(texture, pixels);

                Rectangle source = {0, 0, (float)nx, (float)ny};
                Rectangle dest   = {
                    0, 0, (float)screenWidth, (float)screenHeight
                };
                DrawTexturePro(texture, source, dest, {0, 0}, 0.0f, WHITE);
            } catch (const std::bad_cast &e) {
                std::cout << "Cast failed: " << e.what() << std::endl;
            }
        }

        if (showGrid) {
            Color gridColor = Fade(DARKGRAY, 0.8f); // 使用半透明的深灰色

            auto lines = sim.getLines();
            for (auto line : lines) {
                DrawLine(
                    line.x1 * scaleX,
                    line.y1 * scaleY,
                    line.x2 * scaleX,
                    line.y2 * scaleY,
                    gridColor
                );
            }
        }

        if (showParticle) {
            auto &particles = sim.getParticles();
            for (auto &p : particles) {
                float px = p.x * scaleX;
                float py = p.y * scaleY;

                DrawPixel((int)px, (int)py, RED);
            }
        }
    }
};
