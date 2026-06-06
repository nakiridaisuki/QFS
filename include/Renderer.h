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

    bool showGrid     = true;
    bool showPhi      = false;
    bool showParticle = false;
    bool showVelocity = false;

    void
    DrawVectorArrow(Vector2 start, Vector2 velocity, float scale, Color color) {
        float speedSq = velocity.x * velocity.x + velocity.y * velocity.y;
        if (speedSq < 1e-4f)
            return; // 忽略接近 0 的極小速度

        // 計算箭頭終點
        Vector2 end = {
            start.x + velocity.x * scale, start.y + velocity.y * scale
        };

        // 1. 繪製箭頭主幹
        DrawLineEx(start, end, 1.f, color);
    }

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

    bool &getShowPhi() { return showPhi; }
    bool &getShowGrid() { return showGrid; }
    bool &getShowParticle() { return showParticle; }
    bool &getShowVelocity() { return showVelocity; }

    void draw() {
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

        try {
            const QTSimulator &qt_sim = dynamic_cast<const QTSimulator &>(sim);

            if (showPhi) {
                const float gamma   = 2.2f;
                const float max_val = (float)nx / 2.f;

                for (int j = 0; j < ny; j++) {
                    for (int i = 0; i < nx; i++) {
                        int idx = i + j * nx;
                        float phi =
                            qt_sim.getNodeAt((float)i + 0.5f, (float)j + 0.5f)
                                .phi;

                        Color color;
                        float fraction =
                            1.f - std::clamp(std::abs(phi) / max_val, 0.f, 1.f);
                        float gamma_fraction = std::pow(fraction, gamma);

                        unsigned char value =
                            std::clamp(int(fraction * 255.f), 0, 255);
                        if (phi <= 0.f)
                            color = Color{0, value, 0, 255};
                        else
                            color = Color{value, 0, 0, 255};

                        pixels[idx] = color;
                    }
                }
                UpdateTexture(texture, pixels);

                Rectangle source = {0, 0, (float)nx, (float)ny};
                Rectangle dest   = {
                    0, 0, (float)screenWidth, (float)screenHeight
                };
                DrawTexturePro(texture, source, dest, {0, 0}, 0.0f, WHITE);
            }

            if (showVelocity) {
                auto &ufield     = qt_sim.getUVs(true);
                float velo_scale = 0.1f;
                for (auto &e : ufield) {
                    float px = e.x * scaleX;
                    float py = e.y * scaleY;

                    DrawVectorArrow({px, py}, {e.val, 0}, velo_scale, GREEN);
                }

                auto &vfield = qt_sim.getUVs(false);
                for (auto &e : vfield) {
                    float px = e.x * scaleX;
                    float py = e.y * scaleY;

                    DrawVectorArrow({px, py}, {0, e.val}, velo_scale, GREEN);
                }
            }
        } catch (const std::bad_cast &e) {
            std::cout << "Cast failed: " << e.what() << std::endl;
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
