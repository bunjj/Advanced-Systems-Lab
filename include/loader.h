#include <nlohmann/json.hpp>

#include "geometry.h"

vec load_pos(nlohmann::json& j);
vec load_rot(nlohmann::json& j);
vec load_vec(nlohmann::json& j);
