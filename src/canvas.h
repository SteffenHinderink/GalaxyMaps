#pragma once

#define GL_SILENCE_DEPRECATION

#include "mesh.h"

#include <GLFW/glfw3.h>
#include <nanogui/nanogui.h>

#include <chrono>
#include <fstream>
#include <mutex>

constexpr int SIZE = 888;

enum class Element { CELL,
                     FACE,
                     EDGE,
                     VERTEX,
                     OTHER };

class Canvas : public nanogui::Canvas {
public:
    Canvas(std::string dir, nanogui::Widget* parent, Mesh& mesh, bool parameter_space);

    void set_other(Canvas* other);

    void init();

    void update();

    virtual void draw_contents() override;

    virtual bool mouse_button_event(const nanogui::Vector2i& p, int button, bool down, int modifiers) override;

    virtual bool mouse_drag_event(const nanogui::Vector2i& p, const nanogui::Vector2i& rel, int button, int modifiers) override;

    virtual bool mouse_motion_event(const nanogui::Vector2i& p, const nanogui::Vector2i& rel, int button, int modifiers) override;

    virtual bool scroll_event(const nanogui::Vector2i& p, const nanogui::Vector2f& rel) override;

    virtual bool keyboard_event(int key, int scancode, int action, int modifiers) override;

    std::mutex& get_mutex();

private:
    static nanogui::Matrix4f to_nanogui(Eigen::Matrix4d m);

    Mesh& mesh_;

    bool parameter_space_;

    Canvas* other_;

    std::chrono::steady_clock::time_point last_frame_;

    nanogui::Shader* shader_;

    double scale_;

    Eigen::Matrix4d transformation_matrix_;

    std::vector<int> indices_;

    std::vector<float> positions_;

    std::vector<float> normals_;

    std::vector<float> colors_;

    std::vector<std::pair<Element, int>> sources_;

    bool rot_;

    bool rot_fps_;

    bool trans_;

    int x0_;

    int y0_;

    double size_;

    bool show_cells_;

    bool show_faces_;

    bool show_edges_;

    bool show_vertices_;

    bool show_boundary_;

    bool show_coords_;

    bool show_colored_;

    bool white_background_;

    bool select_;

    OpenVolumeMesh::CellPropertyT<bool> selected_c_;

    OpenVolumeMesh::FacePropertyT<bool> selected_f_;

    OpenVolumeMesh::EdgePropertyT<bool> selected_e_;

    OpenVolumeMesh::VertexPropertyT<bool> selected_v_;

    std::mutex mutex_;

    std::mutex update_mutex_;
};
