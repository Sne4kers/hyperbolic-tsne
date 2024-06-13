#pragma once
#include "cell.hpp"
#include "point.hpp"
#include "math_utils.hpp"
#include <stddef.h>
#include <vector>
#include <algorithm>

class InfinityQuadTree {
public: 
    InfinityQuadTree(std::vector<Point>& points) {
        //std::cout.precision(12);
        //std::cout << "BEGIN BUILDING" << std::endl;
        rec_build_tree(points.begin(), points.end(), Point{-1, -1}, Point{1, 1}, 0);
        //std::cout << "END BUILDING" << std::endl;
    }
    InfinityQuadTree() {}

    std::vector<Cell> get_nodes() {
        return _nodes;
    }

    size_t approximate_centers_of_mass(double x, double y, double theta_sq, double* combined_results) const {
        return approximate_centers_of_mass(Point{x,y}, _nodes.size() - 1, theta_sq, combined_results, 0);
        // std::cout << combined_results.size() << std::endl;
    }

private:
    int rec_build_tree(std::vector<Point>::iterator begin_points, std::vector<Point>::iterator end_points, const Point& min_bounds, const Point& max_bounds, size_t depth){
        if (begin_points == end_points)
            return -1;

        Point bb_center{(max_bounds.x + min_bounds.x) / 2, (max_bounds.y + min_bounds.y) / 2};

        size_t to_be_used = end_points - begin_points;
        // If only 1 point left - it's a leaf
        if (begin_points + 1 == end_points) {
            // Unless the square is not within a unit circle completely
            if (isBoxWithinUnitCircle(min_bounds, max_bounds)) {
                size_t result_idx = _nodes.size();
                _nodes.emplace_back(Cell(depth, min_bounds, max_bounds));

                _nodes[result_idx].is_leaf = true;
                _nodes[result_idx].barycenter = Point{(*begin_points).x, (*begin_points).y};
                _nodes[result_idx].lorentz_factor = hyperbolic_utils::lorentz_factor(_nodes[result_idx].barycenter.to_klein().sq_norm()); 
                _nodes[result_idx].cumulative_size = 1;
                return result_idx;
            }
        }

        // Split the points based on their location
        auto split_y = std::partition(begin_points, end_points, [bb_center](Point a){ return a.y < bb_center.y; });
        auto split_x_lower = std::partition(begin_points, split_y, [bb_center](Point a){ return a.x < bb_center.x; });
        auto split_x_upper = std::partition(split_y, end_points, [bb_center](Point a){ return a.x < bb_center.x; });

        // Recursively call on created partitions
        int child0_idx = rec_build_tree(split_y, split_x_upper, Point{min_bounds.x, bb_center.y}, Point{bb_center.x, max_bounds.y}, depth + 1);
        int child1_idx = rec_build_tree(split_x_upper, end_points, bb_center, max_bounds, depth + 1);
        int child2_idx = rec_build_tree(begin_points, split_x_lower, min_bounds, bb_center, depth + 1);
        int child3_idx = rec_build_tree(split_x_lower, split_y, Point{bb_center.x, min_bounds.y}, Point{max_bounds.x, bb_center.y}, depth + 1);
    
        int only_child = std::max(std::max(child0_idx, child1_idx), std::max(child2_idx, child3_idx));

        if ((child0_idx + child1_idx + child2_idx + child3_idx) == (only_child - 3)) {
            return only_child;
        }

        int result_idx = _nodes.size();
        _nodes.emplace_back(Cell(depth, min_bounds, max_bounds));

        _nodes[result_idx].children_idx[0] = child0_idx;
        _nodes[result_idx].children_idx[1] = child1_idx;
        _nodes[result_idx].children_idx[2] = child2_idx;
        _nodes[result_idx].children_idx[3] = child3_idx;

        // If child_idx is 0, it means that there is no child in that sector of the cell
        double new_lorentz_factor = (child0_idx == -1 ? 0 : _nodes[child0_idx].lorentz_factor)
            + (child1_idx == -1 ? 0 : _nodes[child1_idx].lorentz_factor)
            + (child2_idx == -1 ? 0 : _nodes[child2_idx].lorentz_factor)
            + (child3_idx == -1 ? 0 : _nodes[child3_idx].lorentz_factor);

        Point new_barycenter_klein = ((child0_idx == -1 ? Point{0, 0} : (_nodes[child0_idx].barycenter.to_klein() * _nodes[child0_idx].lorentz_factor))
            + (child1_idx == -1 ? Point{0, 0} : (_nodes[child1_idx].barycenter.to_klein() * _nodes[child1_idx].lorentz_factor))
            + (child2_idx == -1 ? Point{0, 0} : (_nodes[child2_idx].barycenter.to_klein() * _nodes[child2_idx].lorentz_factor))
            + (child3_idx == -1 ? Point{0, 0} : (_nodes[child3_idx].barycenter.to_klein() * _nodes[child3_idx].lorentz_factor))) / new_lorentz_factor;

        double new_max_distance_within_cell = 0;

        for(int i = 0; i < 4; ++i) {
            for(int j = i + 1; j < 4; ++j) {
                if (_nodes[result_idx].children_idx[i] == -1 || _nodes[result_idx].children_idx[j] == -1) {
                    continue;
                }
                new_max_distance_within_cell = std::max(new_max_distance_within_cell, 
                _nodes[_nodes[result_idx].children_idx[i]].barycenter.distance_to_point_poincare(_nodes[_nodes[result_idx].children_idx[j]].barycenter)
                );
            }
        }
        _nodes[result_idx].max_distance_within_squared = new_max_distance_within_cell * new_max_distance_within_cell;

        _nodes[result_idx].barycenter = new_barycenter_klein.to_poincare();
        _nodes[result_idx].lorentz_factor = new_lorentz_factor;
        _nodes[result_idx].cumulative_size = to_be_used;
        _nodes[result_idx].is_leaf = false;
        return result_idx;
    }

    static bool isBoxWithinUnitCircle(const Point& min_bounds, const Point& max_bounds) {
        return hyperbolic_utils::isBoxWithinUnitCircle(min_bounds.x, min_bounds.y, max_bounds.x, max_bounds.y);
    }

    size_t approximate_centers_of_mass(const Point& target, const int& cell_idx, const double& theta_sq, double* combined_results, size_t idx) const {
        auto& current_cell = _nodes[cell_idx];
        //std::cout << "cell_idx " << cell_idx << " idx " << idx << std::endl;

        if (current_cell.is_leaf && std::fabs(target.x - current_cell.barycenter.x) < 1e-5 && std::fabs(target.y - current_cell.barycenter.y) < 1e-5) 
            return idx;
        

        //std::cout << "HERE 1";
        double distance_to_target = target.distance_to_point_poincare(current_cell.barycenter);
        double distance_squared = distance_to_target * distance_to_target;

        //std::cout << "GOT TO HERE" << std::endl;
        // Check the stop condition
        if (current_cell.is_leaf || (!current_cell.contains_infinity && ((current_cell.max_distance_within_squared / distance_squared) < theta_sq))) {
            hyperbolic_utils::distance_grad(target.x, target.y, current_cell.barycenter.x, current_cell.barycenter.y, combined_results[idx], combined_results[idx + 1]);
            combined_results[idx + 2] = distance_to_target;
            combined_results[idx + 3] = current_cell.cumulative_size;
            return idx + 4;
        }
        //std::cout << "LEFT" << std::endl;
        // If stop condition wasn't triggered - go deeper and combine results
        for(int i = 0; i < 4; ++i) {
            if (current_cell.children_idx[i] != -1) {
                //std::cout << "CELL HAS " << current_cell.children_idx[0] << " " << current_cell.children_idx[1] << " " << current_cell.children_idx[2] << " " << current_cell.children_idx[3] << std::endl;
                idx = approximate_centers_of_mass(target, current_cell.children_idx[i], theta_sq, combined_results, idx);
            }
        }
        return idx;
    }

    std::vector<Cell> _nodes;
};