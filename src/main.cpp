#include "Base.h"
#include "GUI.h"
#include "MAC.h"
#include "QT.h"
#include "Renderer.h"
#include <algorithm>
#include <cmath>
#include <vector>

#define RAYGUI_IMPLEMENTATION
#include "raygui.h"

int main() {
    const int screenWidth  = 1000;
    const int screenHeight = screenWidth;
    const int FPS          = 60;

    const std::vector<int> resolutions = {64, 128, 256, 512, 1024};

    int simWidth  = 512;
    int simHeight = simWidth;

    InitWindow(
        screenWidth, screenHeight, "MAC Grid Fluid Simulation - Stable Fluids"
    );
    SetTargetFPS(FPS);

    int simType = SimType::MAC;

    // BaseSimulator *sim = new MACSimulator(simWidth, simHeight);
    BaseSimulator *sim = new QTSimulator(simWidth, simHeight);
    FluidRenderer renderer(sim, screenWidth, screenHeight);
    SimulatorUI gui(sim, renderer, resolutions);

    // 2. 設定 raygui 的全域字體大小與樣式 (選擇性)
    GuiSetStyle(DEFAULT, TEXT_SIZE, 16);

    sim->addWater(sim->getWidth() / 2.0f, sim->getHeight() * 0.5f, 5.0f);

    // === UI 需要的參數變數 ===
    float &brushRadius = gui.brushRadius;
    float &speed       = gui.speed;
    // 定義一塊 UI 區域，用來防止「點擊 UI 時不小心畫出流體」

    int frame_cnt = 0;
    while (!WindowShouldClose()) {
        frame_cnt++;
        Vector2 mousePos = GetMousePosition();

        // 判斷滑鼠是不是在 UI 面板上
        if (frame_cnt % (int)speed == 0) {
            bool isMouseOnUI = gui.inUIPanel(mousePos);

            // 3. 處理滑鼠輸入 (只有滑鼠不在 UI 上時才加水)
            if (IsMouseButtonDown(MOUSE_BUTTON_LEFT) && !isMouseOnUI) {
                float gridX =
                    mousePos.x / (screenWidth / (float)sim->getWidth());
                float gridY =
                    mousePos.y / (screenHeight / (float)sim->getHeight());

                // 這裡原本寫死的 4.0f 改成 UI 變數 brushRadius
                sim->addWater(gridX, gridY, brushRadius);
            }
            if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT) && !isMouseOnUI) {
                float gridX =
                    mousePos.x / (screenWidth / (float)sim->getWidth());
                float gridY =
                    mousePos.y / (screenHeight / (float)sim->getHeight());

                // 這裡原本寫死的 4.0f 改成 UI 變數 brushRadius
                sim->delWater(gridX, gridY, brushRadius);
            }

            float frameTime = std::min(GetFrameTime(), 0.0333f);
            sim->update(frameTime);
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
        gui.drawUI();

        if (gui.need_reset) {
            int new_size = gui.new_sim_size;
            int new_type = gui.new_sim_type;
            if (new_size == simWidth && new_type == simType) {
                sim->reset();
            } else {

                simType  = new_type;
                simWidth = simHeight = resolutions[new_size];

                delete sim;
                if (simType == SimType::MAC)
                    sim = new MACSimulator(simWidth, simHeight);
                else
                    sim = new QTSimulator(simWidth, simHeight);

                renderer.setSimulator(sim);
                gui.setSimulator(sim);
            }
            gui.need_reset = false;
        }

        EndDrawing();
    }

    CloseWindow();
    return 0;
}
