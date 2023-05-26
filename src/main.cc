#include "galaxy.h"
#include "util.h"

#ifdef GUI
#include "canvas.h"

#include <mutex>
#include <thread>

class Screen : public nanogui::Screen {
public:
    Screen() : nanogui::Screen({2 * SIZE, SIZE}, "Galaxy Maps", false) {}

    virtual bool keyboard_event(int key, int scancode, int action, int modifiers) override {
        return callback_(key, scancode, action, modifiers);
    }

    void set_callback(std::function<bool(int, int, int, int)> callback) {
        callback_ = callback;
    }

private:
    std::function<bool(int, int, int, int)> callback_;
};
#endif

void print_usage(std::string dir) {
    std::cout << "Usage: " << dir << "galaxy_map <mesh> <parameter_mesh> [-t] [-o <output>]" << std::endl;
    std::cout << "  e.g. " << dir << "galaxy_map " << dir << "../meshes/example_object.vtk " << dir << "../meshes/example_parameter.vtk" << std::endl;
}

int main(int argc, char** argv) {
    std::string arg0(argv[0]);
    std::string dir = "";
    int pos = std::string(arg0).rfind('/');
    if (pos != std::string::npos) {
        dir = arg0.substr(0, pos + 1);
    }
    std::string input_mesh_name = "";
    std::string parameter_mesh_name = "";
    bool use_tutte = false;
    std::string output_name = "";
    for (int i = 1; i < argc; i++) {
        std::string arg(argv[i]);
        if (arg.length() == 2 && arg.at(0) == '-') {
            char option = arg.at(1);
            switch (option) {
            case 't':
                use_tutte = true;
                break;
            case 'o':
                output_name = std::string(argv[++i]);
                break;
            default:
                print_usage(dir);
                return 0;
            }
        } else if (input_mesh_name == "") {
            input_mesh_name = arg;
        } else if (parameter_mesh_name == "") {
            parameter_mesh_name = arg;
        }
    }

    if (input_mesh_name == "" || parameter_mesh_name == "") {
        print_usage(dir);
        return 0;
    }
    int orientation = 0;
    Mesh mesh = open(input_mesh_name, orientation);
    Mesh parameter_mesh = open(parameter_mesh_name, orientation);
    if (!parameter_from_mesh(mesh, parameter_mesh)) {
        std::cout << "Incompatible input and parameter mesh" << std::endl;
        return 0;
    }
    if (use_tutte) {
        if (!tutte3d(mesh)) {
            std::cout << "3D Tutte failed" << std::endl;
            return 0;
        }
    }

#ifdef GUI
    nanogui::init();
    Screen* screen = new Screen();
    nanogui::Widget* canvases = new nanogui::Widget(screen);
    canvases->set_layout(new nanogui::BoxLayout(nanogui::Orientation::Horizontal, nanogui::Alignment::Fill, 0, 0));
    Canvas* object_canvas = new Canvas(dir, canvases, mesh, false);
    Canvas* parameter_canvas = new Canvas(dir, canvases, mesh, true);
    object_canvas->set_other(parameter_canvas);
    parameter_canvas->set_other(object_canvas);
    std::mutex& object_mutex = object_canvas->get_mutex();
    std::mutex& parameter_mutex = parameter_canvas->get_mutex();
    auto callback = [object_canvas, parameter_canvas](int key, int scancode, int action, int modifiers) {
        bool a = object_canvas->keyboard_event(key, scancode, action, modifiers);
        bool b = parameter_canvas->keyboard_event(key, scancode, action, modifiers);
        return a || b;
    };
    screen->set_callback(callback);
    screen->perform_layout();
    screen->set_visible(true);
    screen->draw_all();
    object_mutex.lock();
    parameter_mutex.lock();
    std::function<void()> run = [&]() {
        galaxy_map(mesh);
        object_mutex.unlock();
        parameter_mutex.unlock();
        object_canvas->update();
        parameter_canvas->update();
    };
    std::thread thread(run);
    nanogui::mainloop(10);
    nanogui::shutdown();
    thread.join();
#else
    galaxy_map(mesh);
#endif

    if (output_name != "") {
        save(mesh, output_name + "_object.vtk", output_name + "_parameter.vtk");
    }
    return 0;
}
