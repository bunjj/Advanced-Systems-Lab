#include "loader.h"
using json = nlohmann::json;

vec load_pos(json& j) {
    return load_vec(j["position"]);
}

vec load_rot(json& j) {
    return load_vec(j["rotation"]);
}

vec load_vec(json& j) { 
    return {j["x"], j["y"], j["z"]}; 
}
