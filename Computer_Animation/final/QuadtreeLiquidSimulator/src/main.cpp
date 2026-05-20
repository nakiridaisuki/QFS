#include "MAC.h"
#include "raylib.h"
#include <algorithm>
#include <iostream>

// 將 Renderer 直接寫在 main 檔案中，保持專案輕量
class FluidRenderer {
  private:
    const MACSimulator &sim;
    int screenWidth, screenHeight;
    float cellWidth, cellHeight;

  public:
    FluidRenderer(const MACSimulator &sim, int screenW, int screenH)
        : sim(sim), screenWidth(screenW), screenHeight(screenH) {
        cellWidth = (float)screenWidth / sim.getWidth();
        cellHeight = (float)screenHeight / sim.getHeight();
    }

    void draw() {
        const auto &dye = sim.getDye();
        int nx = sim.getWidth();
        int ny = sim.getHeight();

        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                // 讀取這格的染料濃度
                double d = dye[i + j * nx];
                if (d > 0.01) {
                    // 將濃度映射到顏色 (我們畫出科技感的藍綠色螢光水流)
                    unsigned char intensity =
                        (unsigned char)std::min(255.0, d * 50.0);
                    Color c = {
                        0,
                        intensity,
                        (unsigned char)std::min(255, intensity + 50),
                        255
                    };

                    DrawRectangle(
                        (int)(i * cellWidth),
                        (int)(j * cellHeight),
                        (int)cellWidth + 1,
                        (int)cellHeight + 1,
                        c
                    );
                }
            }
        }
    }
};

int main() {
    const int screenWidth = 800;
    const int screenHeight = 600;

    // 建立 100x75 的物理網格 (數字越大算越慢，但畫面越細緻)
    MACSimulator sim(120, 90);
    FluidRenderer renderer(sim, screenWidth, screenHeight);

    InitWindow(
        screenWidth, screenHeight, "MAC Grid Fluid Simulation - Stable Fluids"
    );
    SetTargetFPS(60);

    Vector2 prevMousePos = GetMousePosition();

    while (!WindowShouldClose()) {
        float dt = GetFrameTime();

        // 1. 處理滑鼠輸入
        Vector2 mousePos = GetMousePosition();
        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
            // 計算滑鼠移動的向量 (力道)
            float dx = mousePos.x - prevMousePos.x;
            float dy = mousePos.y - prevMousePos.y;

            // 將螢幕座標轉換為「物理網格的索引座標」
            float gridX = mousePos.x / (screenWidth / (float)sim.getWidth());
            float gridY = mousePos.y / (screenHeight / (float)sim.getHeight());

            sim.addForce(gridX, gridY, dx, dy, 4.0f); // 半徑 4.0 格
        }
        prevMousePos = mousePos;

        // 2. 更新物理引擎
        sim.update(dt);

        // 3. 渲染畫面
        BeginDrawing();
        ClearBackground(BLACK); // 黑底讓螢光染料更明顯

        renderer.draw();

        DrawText(
            "Click and Drag to interact with fluid!", 10, 10, 20, LIGHTGRAY
        );
        DrawFPS(10, 35);
        EndDrawing();
    }

    CloseWindow();
    return 0;
}
