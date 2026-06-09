#pragma once

#include <cstdio> // For FILE, popen, pclose
#include <string>
#ifdef _WIN32
#define popen _popen
#define pclose _pclose
#endif

#define RAYGUI_IMPLEMENTATION
#include "raygui.h"

class Recorder {
  private:
    // === 錄影相關全域/靜態變數 ===
    FILE *ffmpeg_pipe = nullptr;
    bool is_recording = false;

    int fps;

  public:
    Recorder(int fps) : fps(fps) {}

    ~Recorder() {
        if (is_recording)
            StopRecording();
    }

    bool isRecording() { return is_recording; }

    void StartRecording(const char *filename) {
        int render_width  = GetRenderWidth();
        int render_height = GetRenderHeight();

        // FFmpeg 指令參數：
        // -y (覆寫檔案), -f rawvideo (輸入格式為原始像素資料), -pix_fmt rgba
        // (顏色格式) -s 寬x高 (輸入解析度), -r 幀率 (輸入幀率), -i - (從 stdin
        // 讀取) -c:v libx264 -pix_fmt yuv420p (輸出相容性佳的 H.264 MP4 格式)
        // -vf vflip (垂直翻轉畫面，修正 OpenGL 座標系顛倒問題)
        std::string cmd = "ffmpeg -y -f rawvideo -pix_fmt rgba -s " +
                          std::to_string(render_width) +
                          "x" +
                          std::to_string(render_height) +
                          " -r " +
                          std::to_string(fps) +
                          " -i - -c:v libx264 -pix_fmt yuv420p -crf 21 " +
                          filename;

#ifdef _WIN32
        ffmpeg_pipe = _popen(cmd.c_str(), "wb"); // Windows 下使用二進位寫入
#else
        ffmpeg_pipe = popen(cmd.c_str(), "w");
#endif

        if (ffmpeg_pipe) {
            is_recording = true;
            TraceLog(LOG_INFO, "Recording started: %s", filename);
        } else {
            TraceLog(
                LOG_WARNING,
                "Failed to start FFmpeg recording pipe. Make sure ffmpeg is "
                "installed and in PATH."
            );
        }
    }

    void RecordFrame() {
        if (!is_recording || !ffmpeg_pipe)
            return;

        Image screen = LoadImageFromScreen();

        fwrite(screen.data, 1, screen.width * screen.height * 4, ffmpeg_pipe);

        UnloadImage(screen);
    }

    void StopRecording() {
        if (ffmpeg_pipe) {
#ifdef _WIN32
            _pclose(ffmpeg_pipe);
#else
            pclose(ffmpeg_pipe);
#endif
            ffmpeg_pipe  = nullptr;
            is_recording = false;
            TraceLog(LOG_INFO, "Recording stopped and video saved.");
        }
    }
};
