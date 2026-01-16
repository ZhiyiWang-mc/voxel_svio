#include "mapManagement.h"
#include "utility.h"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <limits>

bool mapManagement::addPointToVoxel(voxelHashMap &voxel_map, std::shared_ptr<mapPoint> map_point_ptr, double voxel_size, int max_num_points_in_voxel, 
	double min_distance_points, pcl::PointCloud<pcl::PointXYZI>::Ptr voxels_history)
{
	Eigen::Vector3d point = map_point_ptr->getPointXYZ(false);

	short kx = static_cast<short>(point[0] / voxel_size);
	short ky = static_cast<short>(point[1] / voxel_size);
	short kz = static_cast<short>(point[2] / voxel_size);

	map_point_ptr->voxel_idx[0] = kx;
	map_point_ptr->voxel_idx[1] = ky;
	map_point_ptr->voxel_idx[2] = kz;

	voxelHashMap::iterator search = voxel_map.find(voxel(kx, ky, kz));

	if(search != voxel_map.end())
	{
		auto &voxel_block = (search.value());

		if(!voxel_block.IsFull())
		{
			double sq_dist_min_to_points = 10 * voxel_size * voxel_size;
			for (int i(0); i < voxel_block.NumPoints(); ++i)
			{
				auto &map_point_ptr_ = voxel_block.points[i];
				double sq_dist = (map_point_ptr_->getPointXYZ(false) - point).squaredNorm();

				if (sq_dist < sq_dist_min_to_points)
				{
					sq_dist_min_to_points = sq_dist;
				}
			}
			if(sq_dist_min_to_points > (min_distance_points * min_distance_points))
			{
				voxel_block.AddPoint(map_point_ptr);
				return true;
			}
			else
			{
				return false;
			}
		}
	}
	else
	{
		voxelBlock block(max_num_points_in_voxel);
		block.AddPoint(map_point_ptr);
		voxel_map[voxel(kx, ky, kz)] = std::move(block);

		// display
        pcl::PointXYZI point_temp;

        point_temp.x = kx >= 0 ? (kx + 0.5) * voxel_size : (kx - 0.5) * voxel_size;
        point_temp.y = ky >= 0 ? (ky + 0.5) * voxel_size : (ky - 0.5) * voxel_size;
        point_temp.z = kz >= 0 ? (kz + 0.5) * voxel_size : (kz - 0.5) * voxel_size;
        point_temp.intensity = 1;

        voxels_history->points.push_back(point_temp);
        // display

		return true;
    }

	return false;
}

void mapManagement::changeHostVoxel(voxelHashMap &voxel_map, std::shared_ptr<mapPoint> map_point_ptr, double voxel_size, int max_num_points_in_voxel, 
	double min_distance_points, pcl::PointCloud<pcl::PointXYZI>::Ptr voxels_history)
{
	Eigen::Vector3d point = map_point_ptr->getPointXYZ(false);

	short kx = static_cast<short>(point[0] / voxel_size);
	short ky = static_cast<short>(point[1] / voxel_size);
	short kz = static_cast<short>(point[2] / voxel_size);

	if (kx == map_point_ptr->voxel_idx[0] && ky == map_point_ptr->voxel_idx[1] && kz == map_point_ptr->voxel_idx[2]) return;

	voxelHashMap::iterator search_old = voxel_map.find(voxel(map_point_ptr->voxel_idx[0], map_point_ptr->voxel_idx[1], map_point_ptr->voxel_idx[2]));

	auto &voxel_block_old = (search_old.value());
	auto it_old = voxel_block_old.points.begin();

	bool has_found = false;

	while (it_old != voxel_block_old.points.end())
	{
		if ((*it_old) == map_point_ptr)
		{
			has_found = true;
			it_old = voxel_block_old.points.erase(it_old);
			break;
		}

		it_old++;
	}

	assert(has_found);

	map_point_ptr->voxel_idx[0] = kx;
	map_point_ptr->voxel_idx[1] = ky;
	map_point_ptr->voxel_idx[2] = kz;

	voxelHashMap::iterator search_new = voxel_map.find(voxel(kx, ky, kz));

	if(search_new != voxel_map.end())
	{
		auto &voxel_block_new = (search_new.value());
		voxel_block_new.AddPoint(map_point_ptr);
	}
	else
	{
		voxelBlock block_temp(max_num_points_in_voxel);
		block_temp.AddPoint(map_point_ptr);
		voxel_map[voxel(kx, ky, kz)] = std::move(block_temp);

		// display
        pcl::PointXYZI point_temp;

        point_temp.x = kx >= 0 ? (kx + 0.5) * voxel_size : (kx - 0.5) * voxel_size;
        point_temp.y = ky >= 0 ? (ky + 0.5) * voxel_size : (ky - 0.5) * voxel_size;
        point_temp.z = kz >= 0 ? (kz + 0.5) * voxel_size : (kz - 0.5) * voxel_size;
        point_temp.intensity = 1;

        voxels_history->points.push_back(point_temp);
        // display
	}
}

void mapManagement::deleteFromVoxel(voxelHashMap &voxel_map, std::shared_ptr<mapPoint> map_point_ptr)
{
	voxelHashMap::iterator search = voxel_map.find(voxel(map_point_ptr->voxel_idx[0], map_point_ptr->voxel_idx[1], map_point_ptr->voxel_idx[2]));

	auto &voxel_block = (search.value());
	auto it = voxel_block.points.begin();

	bool has_found = false;

	while (it != voxel_block.points.end())
	{
		if (*it == map_point_ptr)
		{
			has_found = true;
			it = voxel_block.points.erase(it);
			break;
		}

		it++;
	}

	assert(has_found);
}

bool mapManagement::updateVoxelLine(voxelHashMap &voxel_map, const voxel &vox, int radius, int min_points, double min_anisotropy,
    double voxel_size, double timestamp)
{
	voxelHashMap::iterator search = voxel_map.find(voxel(vox.x, vox.y, vox.z));
	if (search == voxel_map.end())
		return false;

	auto &voxel_block = (search.value());

	if (voxel_block.line_update_time >= timestamp)
		return voxel_block.line_quality > 0.0;

	std::vector<Eigen::Vector3d> pts;

	for (int dx = -radius; dx <= radius; dx++)
	{
		for (int dy = -radius; dy <= radius; dy++)
		{
			for (int dz = -radius; dz <= radius; dz++)
			{
				voxelHashMap::iterator search_neighbor = voxel_map.find(voxel(vox.x + dx, vox.y + dy, vox.z + dz));
				if (search_neighbor == voxel_map.end())
					continue;

				auto &block_neighbor = (search_neighbor.value());
				for (const auto &mp : block_neighbor.points)
				{
					pts.push_back(mp->getPointXYZ(false));
				}
			}
		}
	}

	if ((int)pts.size() < min_points)
	{
		voxel_block.line_quality = 0.0;
		voxel_block.line_dir.setZero();
		voxel_block.line_point.setZero();
		voxel_block.line_update_time = timestamp;
		return false;
	}

	double fit_dist = 0.25 * voxel_size;
	if (setting_line_voxel_max_dist > 0.0f)
		fit_dist = std::min(fit_dist, 0.5 * (double)setting_line_voxel_max_dist);
	fit_dist = std::max(fit_dist, 0.02);

	int max_iter = std::min(120, std::max(30, (int)pts.size() * 6));
	int best_inliers = 0;
	Eigen::Vector3d best_point = Eigen::Vector3d::Zero();
	Eigen::Vector3d best_dir = Eigen::Vector3d::Zero();

	int iter = 0;
	for (size_t i = 0; i < pts.size(); i++)
	{
		for (size_t j = i + 1; j < pts.size(); j++)
		{
			if (iter++ >= max_iter)
				goto finish_ransac;

			Eigen::Vector3d dir = pts[j] - pts[i];
			double norm = dir.norm();
			if (norm < 1e-6)
				continue;

			dir /= norm;
			int inliers = 0;
			for (const auto &p : pts)
			{
				double dist = (p - pts[i]).cross(dir).norm();
				if (dist <= fit_dist)
					inliers++;
			}

			if (inliers > best_inliers)
			{
				best_inliers = inliers;
				best_point = pts[i];
				best_dir = dir;
			}
		}
	}

finish_ransac:
	if (best_inliers < min_points)
	{
		voxel_block.line_quality = 0.0;
		voxel_block.line_dir.setZero();
		voxel_block.line_point.setZero();
		voxel_block.line_update_time = timestamp;
		return false;
	}

	std::vector<Eigen::Vector3d> inliers;
	inliers.reserve(best_inliers);
	for (const auto &p : pts)
	{
		double dist = (p - best_point).cross(best_dir).norm();
		if (dist <= fit_dist)
			inliers.push_back(p);
	}

	double inlier_ratio = (double)inliers.size() / std::max(1.0, (double)pts.size());
	if ((int)inliers.size() < min_points || inlier_ratio < 0.6)
	{
		voxel_block.line_quality = 0.0;
		voxel_block.line_dir.setZero();
		voxel_block.line_point.setZero();
		voxel_block.line_update_time = timestamp;
		return false;
	}

	Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
	for (const auto &p : inliers)
		centroid += p;
	centroid /= (double)inliers.size();

	Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
	for (const auto &p : inliers)
	{
		Eigen::Vector3d d = p - centroid;
		cov.noalias() += d * d.transpose();
	}
	cov /= std::max(1, (int)inliers.size() - 1);

	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(cov);
	if (solver.info() != Eigen::Success)
	{
		voxel_block.line_quality = 0.0;
		voxel_block.line_dir.setZero();
		voxel_block.line_point = centroid;
		voxel_block.line_update_time = timestamp;
		return false;
	}

	Eigen::Vector3d evals = solver.eigenvalues();
	Eigen::Matrix3d evecs = solver.eigenvectors();

	Eigen::Index max_idx;
	evals.maxCoeff(&max_idx);

	double lambda_max = evals(max_idx);
	double lambda_second = -1.0;
	for (int i = 0; i < 3; i++)
	{
		if (i == max_idx)
			continue;
		if (evals(i) > lambda_second)
			lambda_second = evals(i);
	}
	double anisotropy = (lambda_max - lambda_second) / (lambda_max + 1e-6);

	if (anisotropy < min_anisotropy)
	{
		voxel_block.line_quality = 0.0;
		voxel_block.line_dir.setZero();
		voxel_block.line_point = centroid;
		voxel_block.line_update_time = timestamp;
		return false;
	}

	Eigen::Vector3d line_dir = evecs.col(max_idx);
	if (line_dir.norm() < 1e-6)
	{
		voxel_block.line_quality = 0.0;
		voxel_block.line_dir.setZero();
		voxel_block.line_point = centroid;
		voxel_block.line_update_time = timestamp;
		return false;
	}

	line_dir.normalize();

	double min_proj = std::numeric_limits<double>::infinity();
	double max_proj = -std::numeric_limits<double>::infinity();
	for (const auto &p : inliers)
	{
		double proj = line_dir.dot(p - centroid);
		if (proj < min_proj)
			min_proj = proj;
		if (proj > max_proj)
			max_proj = proj;
	}

	double span = max_proj - min_proj;
	double min_span = 0.3 * voxel_size;
	if (span < min_span)
	{
		voxel_block.line_quality = 0.0;
		voxel_block.line_dir.setZero();
		voxel_block.line_point = centroid;
		voxel_block.line_update_time = timestamp;
		return false;
	}

	voxel_block.line_point = centroid;
	voxel_block.line_dir = line_dir;
	voxel_block.line_quality = anisotropy * inlier_ratio;
	voxel_block.line_update_time = timestamp;

	return true;
}

bool mapManagement::getVoxelLine(const voxelHashMap &voxel_map, const voxel &vox, Eigen::Vector3d &line_point, Eigen::Vector3d &line_dir, double &line_quality)
{
	voxelHashMap::const_iterator search = voxel_map.find(voxel(vox.x, vox.y, vox.z));
	if (search == voxel_map.end())
		return false;

	const auto &voxel_block = (search.value());
	if (voxel_block.line_quality <= 0.0 || voxel_block.line_dir.squaredNorm() < 1e-9)
		return false;

	line_point = voxel_block.line_point;
	line_dir = voxel_block.line_dir;
	line_quality = voxel_block.line_quality;
	return true;
}
