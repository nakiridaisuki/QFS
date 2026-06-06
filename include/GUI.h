#pragma once

#include "Base.h"
#include "Renderer.h"
#include "raygui.h"
#include <string>
#include <vector>

namespace SimType {
enum SimulatorType { MAC = 0, QT = 1 };
}

class SimulatorUI {
    struct PackedSliderElement {
        const char *textLeft;
        float *value;
        float min_val, max_val;
    };
    struct PackedCheckBoxElement {
        const char *text;
        bool *value;
    };

    float global_y_gap = 2;

  public:
    BaseSimulator *sim;
    FluidRenderer &renderer;
    float brushRadius;
    float gravity;
    float tension;
    bool showUI;
    bool showDebug;
    int iterations;
    float speed;
    int frame_cnt;
    int new_sim_type; // 0: MAC, 1: QT
    int new_sim_size;
    bool need_reset;
    Rectangle uiPanelRec;

    std::vector<PackedSliderElement> main_sliders;
    std::vector<PackedCheckBoxElement> debug_checkbox;

    std::string resolution_data;

    SimulatorUI(
        BaseSimulator *sim,
        FluidRenderer &renderer,
        std::vector<int> resolutions
    )
        : sim(sim), renderer(renderer) {

        // Define init value
        brushRadius = 4.0f;
        gravity     = sim->getGravity();
        tension     = sim->getSurfaceTension();
        showUI      = true;
        showDebug   = true;
        iterations  = 1;
        speed       = 1.f;
        frame_cnt   = 0;
        uiPanelRec  = {10, 60, 260, 150};

        new_sim_size = 1;
        new_sim_type = 0;
        need_reset   = false;

        resolution_data = "";
        for (int i = 0; i < resolutions.size(); i++) {
            if (i != 0)
                resolution_data.push_back(';');
            resolution_data += std::to_string(resolutions[i]);
        }

        // Define main sliders
        main_sliders.push_back({"Brush Size", &brushRadius, 1.0f, 200.0f});
        main_sliders.push_back({"Speed", &speed, 1.0f, 60.0f});
        main_sliders.push_back({"Gravity", &gravity, 0.0f, 1000.0f});
        main_sliders.push_back({"Tension", &tension, 0.0f, 1000.0f});

        // Define debug checkboxes
        debug_checkbox.push_back({"Particles", &renderer.getShowParticle()});
        debug_checkbox.push_back({"Grid", &renderer.getShowGrid()});
        debug_checkbox.push_back({"Phi", &renderer.getShowPhi()});
        debug_checkbox.push_back({"Velocity", &renderer.getShowVelocity()});
    }

    void setSimulator(BaseSimulator *new_sim) { sim = new_sim; }

    bool inUIPanel(Vector2 pos) {
        return showUI && CheckCollisionPointRec(pos, uiPanelRec);
    }

    void drawUI() {
        float curr_y = 60;
        GuiCheckBox(Rectangle{20, curr_y, 20, 20}, "Show UI", &showUI);
        GuiCheckBox(Rectangle{120, curr_y, 20, 20}, "Show Debug", &showDebug);
        curr_y += 20 + global_y_gap;

        if (showUI) {
            // 畫一個半透明的背景板，讓 UI 更清楚
            // DrawRectangleRec(uiPanelRec, Fade(DARKGRAY, 0.8f));

            curr_y = drawPackedSliders(120, curr_y, main_sliders);

            sim->setGravity(gravity);
            sim->setSigma(tension);

            GuiSpinner(
                Rectangle{120, curr_y, 120, 20},
                "Max Iteration",
                &iterations,
                0,
                10,
                false
            );
            updateY(curr_y, 20);

            GuiComboBox(
                Rectangle{20, curr_y, 120, 20},
                resolution_data.c_str(),
                &new_sim_size
            );
            updateY(curr_y, 20);

            GuiComboBox(
                Rectangle{20, curr_y, 120, 20}, "MAC;Quadtree", &new_sim_type
            );
            need_reset =
                GuiButton(Rectangle{150, curr_y, 100, 20}, "Reset Fluid");
            updateY(curr_y, 20);

            if (showDebug)
                curr_y = drawPackedCheckBoxes(20, curr_y, debug_checkbox);
        }
    }

    void updateY(float &y, float h) { y += h + global_y_gap; }

    float drawPackedSliders(
        float baseX, float baseY, std::vector<PackedSliderElement> &elements
    ) {
        const int width  = 120;
        const int height = 20;

        float currY = baseY;

        for (auto &e : elements) {
            GuiSlider(
                Rectangle{baseX, currY, width, height},
                e.textLeft,
                TextFormat("%.1f", *e.value),
                e.value,
                e.min_val,
                e.max_val
            );
            currY += global_y_gap + height;
        }
        return currY;
    }

    float drawPackedCheckBoxes(
        float baseX, float baseY, std::vector<PackedCheckBoxElement> &elements
    ) {
        const int width  = 20;
        const int height = 20;
        const int x_gap  = 100;

        float currY  = baseY;
        bool is_left = true;

        for (auto &e : elements) {
            if (is_left) {
                GuiCheckBox(
                    Rectangle{baseX, currY, width, height}, e.text, e.value
                );
                is_left = false;
            } else {
                GuiCheckBox(
                    Rectangle{baseX + x_gap, currY, width, height},
                    e.text,
                    e.value
                );
                currY += global_y_gap + height;
                is_left = true;
            }
        }

        if (!is_left)
            currY += global_y_gap + height;
        return currY;
    }
};
