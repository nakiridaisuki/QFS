#include "MAC.h"
#include "QT.h"
#include "Renderer.h"
#include <algorithm>
#include <cmath>

#define RAYGUI_IMPLEMENTATION
#include "raygui.h"

int main() {
    const int screenWidth  = 1000;
    const int screenHeight = 1000;
    const int simWidth     = 128;
    const int simHeight    = simWidth;
    const int FPS          = 60;

    InitWindow(
        screenWidth, screenHeight, "MAC Grid Fluid Simulation - Stable Fluids"
    );
    SetTargetFPS(FPS);

    // MACSimulator sim(simWidth, simHeight);
    QTSimulator sim(simWidth, simHeight);
    FluidRenderer renderer(sim, screenWidth, screenHeight);

    // 2. 設定 raygui 的全域字體大小與樣式 (選擇性)
    GuiSetStyle(DEFAULT, TEXT_SIZE, 16);

    sim.addWater(sim.getWidth() / 2.0f, sim.getHeight() * 0.5f, 5.0f);

    // === UI 需要的參數變數 ===
    float brushRadius = 4.0f;
    float gravity     = sim.getGravity();
    float tension     = sim.getSurfaceTension();
    bool showUI       = true;
    int iterations    = 1;
    float speed       = 1.f;
    int frame_cnt     = 0;
    // 定義一塊 UI 區域，用來防止「點擊 UI 時不小心畫出流體」
    Rectangle uiPanelRec = {10, 60, 260, 150};

    while (!WindowShouldClose()) {
        frame_cnt++;
        Vector2 mousePos = GetMousePosition();

        // 判斷滑鼠是不是在 UI 面板上
        if (frame_cnt % (int)speed == 0) {
            bool isMouseOnUI =
                showUI && CheckCollisionPointRec(mousePos, uiPanelRec);

            // 3. 處理滑鼠輸入 (只有滑鼠不在 UI 上時才加水)
            if (IsMouseButtonDown(MOUSE_BUTTON_LEFT) && !isMouseOnUI) {
                float gridX =
                    mousePos.x / (screenWidth / (float)sim.getWidth());
                float gridY =
                    mousePos.y / (screenHeight / (float)sim.getHeight());

                // 這裡原本寫死的 4.0f 改成 UI 變數 brushRadius
                sim.addWater(gridX, gridY, brushRadius);
            }
            if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT) && !isMouseOnUI) {
                float gridX =
                    mousePos.x / (screenWidth / (float)sim.getWidth());
                float gridY =
                    mousePos.y / (screenHeight / (float)sim.getHeight());

                // 這裡原本寫死的 4.0f 改成 UI 變數 brushRadius
                sim.delWater(gridX, gridY, brushRadius);
            }

            float frameTime = std::min(GetFrameTime(), 0.0333f);
            float idel_dt;
            if (sim.getMaxVel() == 0)
                idel_dt = frameTime;
            else
                idel_dt = 1.f / sim.getMaxVel();
            int idel_iter = std::ceil(frameTime / idel_dt);

            int actual_iter = std::clamp(idel_iter, 1, iterations);
            float actual_dt = frameTime / actual_iter;
            for (int iter = 0; iter < actual_iter; iter++) {
                sim.update(actual_dt);
            }
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
                Rectangle{100, 90, 120, 20},
                "Speed",
                TextFormat("%.1f", speed),
                &speed,
                1.0f,
                60.0f
            );

            GuiSlider(
                Rectangle{100, 110, 120, 20},
                "Gravity",
                TextFormat("%.0f", gravity),
                &gravity,
                0.0f,
                1000.0f
            );

            GuiSlider(
                Rectangle{100, 130, 120, 20},
                "Tension",
                TextFormat("%.0f", tension),
                &tension,
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
            sim.setSigma(tension);

            GuiCheckBox(
                Rectangle{10, 170, 20, 20},
                "Particles",
                &renderer.getShowParticle()
            );
            GuiCheckBox(
                Rectangle{10, 190, 20, 20}, "Grid", &renderer.getShowGrid()
            );
            GuiCheckBox(
                Rectangle{110, 170, 20, 20}, "Phi", &renderer.getShowPhi()
            );
            GuiCheckBox(
                Rectangle{110, 190, 20, 20},
                "Velocity",
                &renderer.getShowVelocity()
            );

            // 把 Reset 按鈕稍微往下挪
            if (GuiButton(Rectangle{100, 215, 120, 30}, "Reset Fluid")) {
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
