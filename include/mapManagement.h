#pragma once

// c++ include
#include <iostream>
#include <memory>
#include <mutex>
#include <vector>

// lib include
#include <Eigen/Core>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl_ros/point_cloud.h>
#include <pcl_conversions/pcl_conversions.h>

// function include
#include "mapPoint.h"

class mapManagement
{
public:

  static bool addPointToVoxel(voxelHashMap &voxel_map, std::shared_ptr<mapPoint> map_point_ptr, double voxel_size, int max_num_points_in_voxel, 
      double min_distance_points, pcl::PointCloud<pcl::PointXYZI>::Ptr voxels_history);

  static void changeHostVoxel(voxelHashMap &voxel_map, std::shared_ptr<mapPoint> map_point_ptr, double voxel_size, int max_num_points_in_voxel, 
      double min_distance_points, pcl::PointCloud<pcl::PointXYZI>::Ptr voxels_history);

  static void deleteFromVoxel(voxelHashMap &voxel_map, std::shared_ptr<mapPoint> map_point_ptr);

  static bool updateVoxelLine(voxelHashMap &voxel_map, const voxel &vox, int radius, int min_points, double min_anisotropy,
      double voxel_size, double timestamp);

  static bool getVoxelLine(const voxelHashMap &voxel_map, const voxel &vox, Eigen::Vector3d &line_point, Eigen::Vector3d &line_dir, double &line_quality);

private:

  mapManagement();
};
