#include "MAC.h"
#include "raylib.h"
#include <algorithm>

#define RAYGUI_IMPLEMENTATION
#include "raygui.h"

class FluidRenderer {
  private:
    const MACSimulator &sim;
    int screenWidth, screenHeight;
    Image image;
    Texture2D texture;
    Color *pixels;

  public:
    FluidRenderer(const MACSimulator &sim, int screenW, int screenH)
        : sim(sim), screenWidth(screenW), screenHeight(screenH) {

        // 建立與模擬網格大小相同的 Image
        pixels = new Color[sim.getWidth() * sim.getHeight()];
        image = {
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

    void draw() {
        const auto &celltype = sim.getCell();
        int nx = sim.getWidth();
        int ny = sim.getHeight();

        // 在 CPU 端快速填充像素
        for (int i = 0; i < nx * ny; ++i) {
            if (celltype[i]) {
                pixels[i] = Color{0, 50, 100, 255}; // 你的科技螢光藍
            } else {
                pixels[i] = BLANK; // 透明或黑色
            }
        }

        // 一次性將資料送給 GPU
        UpdateTexture(texture, pixels);

        // 放大畫回螢幕上
        Rectangle source = {0, 0, (float)nx, (float)ny};
        Rectangle dest = {0, 0, (float)screenWidth, (float)screenHeight};
        DrawTexturePro(texture, source, dest, {0, 0}, 0.0f, WHITE);
    }
};

int main() {
    const int screenWidth = 1200;
    const int screenHeight = 900;
    const int simWidth = screenWidth;
    const int simHeight = screenHeight;
    const int FPS = 60;

    InitWindow(
        screenWidth, screenHeight, "MAC Grid Fluid Simulation - Stable Fluids"
    );
    SetTargetFPS(FPS);

    MACSimulator sim(simWidth, simHeight, FPS);
    FluidRenderer renderer(sim, screenWidth, screenHeight);

    // 2. 設定 raygui 的全域字體大小與樣式 (選擇性)
    GuiSetStyle(DEFAULT, TEXT_SIZE, 16);

    sim.addWater(sim.getWidth() / 2.0f, sim.getHeight() * 0.5f, 25.0f);

    // === UI 需要的參數變數 ===
    float brushRadius = 4.0f;
    float gravity = sim.getGravity();
    bool showUI = true;
    int iterations = 1;
    // 定義一塊 UI 區域，用來防止「點擊 UI 時不小心畫出流體」
    Rectangle uiPanelRec = {10, 60, 260, 150};

    while (!WindowShouldClose()) {
        Vector2 mousePos = GetMousePosition();

        // 判斷滑鼠是不是在 UI 面板上
        bool isMouseOnUI =
            showUI && CheckCollisionPointRec(mousePos, uiPanelRec);

        // 3. 處理滑鼠輸入 (只有滑鼠不在 UI 上時才加水)
        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT) && !isMouseOnUI) {
            float gridX = mousePos.x / (screenWidth / (float)sim.getWidth());
            float gridY = mousePos.y / (screenHeight / (float)sim.getHeight());

            // 這裡原本寫死的 4.0f 改成 UI 變數 brushRadius
            sim.addWater(gridX, gridY, brushRadius);
        }
        if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT) && !isMouseOnUI) {
            float gridX = mousePos.x / (screenWidth / (float)sim.getWidth());
            float gridY = mousePos.y / (screenHeight / (float)sim.getHeight());

            // 這裡原本寫死的 4.0f 改成 UI 變數 brushRadius
            sim.delWater(gridX, gridY, brushRadius);
        }

        float frameTime = std::min(GetFrameTime(), 0.0333f);
        float dt = frameTime / (float)iterations;
        for (int iter = 0; iter < iterations; iter++) {
            sim.update(dt);
        }

        // 4. 渲染畫面與 UI
        BeginDrawing();
        ClearBackground(BLACK);

        renderer.draw();

        // 畫一點基本的文字
        DrawText(
            "Click and Drag to interact with fluid!", 10, 10, 20, LIGHTGRAY
        );
        DrawFPS(10, 35);

        // === 繪製 UI 面板 ===
        GuiCheckBox(Rectangle{250, 10, 20, 20}, "Show UI", &showUI);

        if (showUI) {
            // 畫一個半透明的背景板，讓 UI 更清楚
            DrawRectangleRec(uiPanelRec, Fade(DARKGRAY, 0.8f));

            // 加入滑桿 Slider (Bounds, 左邊文字, 右邊文字顯示數值, 變數指標,
            // 最小值, 最大值)
            GuiSlider(
                Rectangle{100, 70, 120, 20},
                "Brush Size",
                TextFormat("%.1f", brushRadius),
                &brushRadius,
                1.0f,
                200.0f
            );

            GuiSlider(
                Rectangle{100, 130, 120, 20},
                "Gravity",
                TextFormat("%.0f", gravity),
                &gravity,
                0.0f,
                1000.0f
            );

            GuiSpinner(
                Rectangle{100, 150, 120, 20},
                "Iteration times",
                &iterations,
                0,
                10,
                false
            );

            sim.setGravity(gravity);

            // 你也可以在這邊加一顆按鈕來重置場景
            if (GuiButton(Rectangle{100, 170, 120, 30}, "Reset Fluid")) {
                sim.reset();
                sim.addWater(
                    sim.getWidth() / 2.0f, sim.getHeight() * 0.5f, 25.0f
                );
            }
        }

        EndDrawing();
    }

    CloseWindow();
    return 0;
}
