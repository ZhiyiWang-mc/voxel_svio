#include "updaterHelper.h"
#include "state.h"
#include "quatOps.h"
#include <algorithm>
#include <cmath>

void updaterHelper::getFeatureJacobianRepresentation(std::shared_ptr<state> state_ptr, updaterHelperFeature &feature, Eigen::MatrixXd &H_f,
	std::vector<Eigen::MatrixXd> &H_x, std::vector<std::shared_ptr<baseType>> &x_order)
{
	H_f.resize(3, 3);
	H_f.setIdentity();
	return;
}

void updaterHelper::getFeatureJacobianFull(std::shared_ptr<state> state_ptr, updaterHelperFeature &feature, Eigen::MatrixXd &H_f, 
	Eigen::MatrixXd &H_x, Eigen::VectorXd &res, std::vector<std::shared_ptr<baseType>> &x_order)
{
	int total_meas = 0;
	for (auto const &pair : feature.timestamps)
	{
		total_meas += (int)pair.second.size();
	}

	bool use_lines = setting_line_enable && (total_meas >= setting_line_min_track);

	int total_hx = 0;
	std::unordered_map<std::shared_ptr<baseType>, size_t> map_hx;
	for (auto const &pair : feature.timestamps)
	{
		std::shared_ptr<poseJpl> calibration = state_ptr->calib_cam_imu.at(pair.first);
		std::shared_ptr<vec> distortion = state_ptr->cam_intrinsics.at(pair.first);

		if (state_ptr->options.do_calib_camera_pose)
		{
			map_hx.insert({calibration, total_hx});
			x_order.push_back(calibration);
			total_hx += calibration->getSize();
		}

		if (state_ptr->options.do_calib_camera_intrinsics)
		{
			map_hx.insert({distortion, total_hx});
			x_order.push_back(distortion);
			total_hx += distortion->getSize();
		}

		for (size_t m = 0; m < feature.timestamps[pair.first].size(); m++)
		{
			std::shared_ptr<poseJpl> clone_Ci = state_ptr->clones_imu.at(feature.timestamps[pair.first].at(m));

			if (map_hx.find(clone_Ci) == map_hx.end())
			{
				map_hx.insert({clone_Ci, total_hx});
				x_order.push_back(clone_Ci);
				total_hx += clone_Ci->getSize();
			}
		}
	}

	Eigen::Vector3d p_FinG = feature.position_global;
	Eigen::Vector3d p_FinG_fej = feature.position_global_fej;

	int jacobsize = 3;
	int total_rows = 2 * total_meas;

	res = Eigen::VectorXd::Zero(total_rows);
	H_f = Eigen::MatrixXd::Zero(total_rows, jacobsize);
	H_x = Eigen::MatrixXd::Zero(total_rows, total_hx);

	Eigen::MatrixXd dpfg_dlambda;
	std::vector<Eigen::MatrixXd> dpfg_dx;
	std::vector<std::shared_ptr<baseType>> dpfg_dx_order;
	updaterHelper::getFeatureJacobianRepresentation(state_ptr, feature, dpfg_dlambda, dpfg_dx, dpfg_dx_order);

	int row = 0;
	for (auto const &pair : feature.timestamps)
	{
		std::shared_ptr<vec> distortion = state_ptr->cam_intrinsics.at(pair.first);
		std::shared_ptr<poseJpl> calibration = state_ptr->calib_cam_imu.at(pair.first);
		Eigen::Matrix3d R_ItoC = calibration->getRot();
		Eigen::Vector3d p_IinC = calibration->getPos();

		for (size_t m = 0; m < feature.timestamps[pair.first].size(); m++)
		{
			std::shared_ptr<poseJpl> clone_Ii = state_ptr->clones_imu.at(feature.timestamps[pair.first].at(m));
			Eigen::Matrix3d R_GtoIi = clone_Ii->getRot();
			Eigen::Vector3d p_IiinG = clone_Ii->getPos();

			Eigen::Vector3d p_FinIi = R_GtoIi * (p_FinG - p_IiinG);

			Eigen::Vector3d p_FinCi = R_ItoC * p_FinIi + p_IinC;
			Eigen::Vector2d uv_norm;
			uv_norm << p_FinCi(0) / p_FinCi(2), p_FinCi(1) / p_FinCi(2);

			Eigen::Vector2d uv_dist;
			uv_dist = state_ptr->cam_intrinsics_cameras.at(pair.first)->distortD(uv_norm);

			Eigen::Vector2d uv_m;
			uv_m << (double)feature.uvs[pair.first].at(m)(0), (double)feature.uvs[pair.first].at(m)(1);
			Eigen::Vector2d r_uv = uv_m - uv_dist;
			res.block(row, 0, 2, 1) = r_uv;

			if (state_ptr->options.do_fej)
			{
				R_GtoIi = clone_Ii->getRotFej();
				p_IiinG = clone_Ii->getPosFej();

				p_FinIi = R_GtoIi * (p_FinG_fej - p_IiinG);
				p_FinCi = R_ItoC * p_FinIi + p_IinC;
			}

			Eigen::MatrixXd dz_dzn, dz_dzeta;
			state_ptr->cam_intrinsics_cameras.at(pair.first)->computeDistortJacobian(uv_norm, dz_dzn, dz_dzeta);

			Eigen::MatrixXd dzn_dpfc = Eigen::MatrixXd::Zero(2, 3);
			dzn_dpfc << 1 / p_FinCi(2), 0, -p_FinCi(0) / (p_FinCi(2) * p_FinCi(2)), 0, 1 / p_FinCi(2), -p_FinCi(1) / (p_FinCi(2) * p_FinCi(2));

			Eigen::MatrixXd dpfc_dpfg = R_ItoC * R_GtoIi;

			Eigen::MatrixXd dpfc_dclone = Eigen::MatrixXd::Zero(3, 6);
			dpfc_dclone.block(0, 0, 3, 3).noalias() = R_ItoC * quatType::skewSymmetric(p_FinIi);
			dpfc_dclone.block(0, 3, 3, 3) = - dpfc_dpfg;

			Eigen::MatrixXd dz_dpfc = dz_dzn * dzn_dpfc;
			Eigen::MatrixXd dz_dpfg = dz_dpfc * dpfc_dpfg;

			Eigen::MatrixXd Hf_point = dz_dpfg * dpfg_dlambda;
			Eigen::MatrixXd Hx_point = dz_dpfc * dpfc_dclone;

			std::vector<Eigen::MatrixXd> Hx_point_extra;
			for (size_t i = 0; i < dpfg_dx_order.size(); i++)
			{
				Hx_point_extra.push_back(dz_dpfg * dpfg_dx.at(i));
			}

			Eigen::MatrixXd dpfc_dcalib;
			if (state_ptr->options.do_calib_camera_pose)
			{
				dpfc_dcalib = Eigen::MatrixXd::Zero(3, 6);
				dpfc_dcalib.block(0, 0, 3, 3) = quatType::skewSymmetric(p_FinCi - p_IinC);
				dpfc_dcalib.block(0, 3, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
			}

			Eigen::MatrixXd Hx_dist;
			if (state_ptr->options.do_calib_camera_intrinsics)
			{
				Hx_dist = dz_dzeta;
			}

			int block_rows = 2;
			bool has_line = false;
			double w_line = 1.0;
			double w_point = 1.0;
			Eigen::RowVector2d nT;
			Eigen::RowVector2d tT;

			if (use_lines)
			{
				auto it_dir = feature.line_directions.find(pair.first);
				auto it_w = feature.line_weights.find(pair.first);

				if (it_dir != feature.line_directions.end() && it_w != feature.line_weights.end() &&
					it_dir->second.size() > m && it_w->second.size() > m)
				{
					Eigen::Vector2d line_dir = it_dir->second.at(m).cast<double>();
					double norm = line_dir.norm();

					if (it_w->second.at(m) > 0.0f && norm > 1e-6)
					{
						line_dir /= norm;
						nT << -line_dir(1), line_dir(0);
						tT << line_dir(0), line_dir(1);
						w_line = sqrt((double)it_w->second.at(m));
						has_line = true;
					}
				}

				if (has_line && setting_line_point_gate > 0.0f)
				{
					double point_res = r_uv.norm();
					if (point_res > setting_line_point_gate)
						has_line = false;
				}

				if (has_line && setting_line_rot_deg > 0.0f && feature.anchor_clone_timestamp >= 0)
				{
					auto it_anchor = state_ptr->clones_imu.find(feature.anchor_clone_timestamp);
					if (it_anchor != state_ptr->clones_imu.end())
					{
						Eigen::Matrix3d R_rel = clone_Ii->getRot() * it_anchor->second->getRot().transpose();
						double cos_theta = 0.5 * (R_rel.trace() - 1.0);
						cos_theta = std::min(1.0, std::max(-1.0, cos_theta));
						double angle_deg = std::acos(cos_theta) * 180.0 / 3.14159265358979323846;

						if (angle_deg > setting_line_rot_deg)
						{
							double line_scale = std::min(1.0, w_line);
							double rot_scale = 1.0 - (1.0 - setting_line_rot_point_weight) * line_scale;
							w_point *= rot_scale;
						}
					}
				}
			}

			if (has_line)
			{
				res(row, 0) = w_point * (tT * r_uv)(0, 0);
				res(row + 1, 0) = w_line * (nT * r_uv)(0, 0);

				Eigen::MatrixXd t_dz_dpfc = tT * dz_dpfc;
				Eigen::MatrixXd n_dz_dpfc = nT * dz_dpfc;
				Eigen::MatrixXd t_dz_dpfg = tT * dz_dpfg;
				Eigen::MatrixXd n_dz_dpfg = nT * dz_dpfg;

				H_f.block(row, 0, 1, H_f.cols()).noalias() = w_point * t_dz_dpfg * dpfg_dlambda;
				H_f.block(row + 1, 0, 1, H_f.cols()).noalias() = w_line * n_dz_dpfg * dpfg_dlambda;
				H_x.block(row, map_hx[clone_Ii], 1, clone_Ii->getSize()).noalias() = w_point * t_dz_dpfc * dpfc_dclone;
				H_x.block(row + 1, map_hx[clone_Ii], 1, clone_Ii->getSize()).noalias() = w_line * n_dz_dpfc * dpfc_dclone;

				for (size_t i = 0; i < dpfg_dx_order.size(); i++)
				{
					H_x.block(row, map_hx[dpfg_dx_order.at(i)], 1, dpfg_dx_order.at(i)->getSize()).noalias() += 
						w_point * t_dz_dpfg * dpfg_dx.at(i);
					H_x.block(row + 1, map_hx[dpfg_dx_order.at(i)], 1, dpfg_dx_order.at(i)->getSize()).noalias() += 
						w_line * n_dz_dpfg * dpfg_dx.at(i);
				}

				if (state_ptr->options.do_calib_camera_pose)
				{
					H_x.block(row, map_hx[calibration], 1, calibration->getSize()).noalias() += w_point * t_dz_dpfc * dpfc_dcalib;
					H_x.block(row + 1, map_hx[calibration], 1, calibration->getSize()).noalias() += w_line * n_dz_dpfc * dpfc_dcalib;
				}

				if (state_ptr->options.do_calib_camera_intrinsics)
				{
					H_x.block(row, map_hx[distortion], 1, distortion->getSize()) = w_point * tT * Hx_dist;
					H_x.block(row + 1, map_hx[distortion], 1, distortion->getSize()) = w_line * nT * Hx_dist;
				}
			}
			else
			{
				H_f.block(row, 0, 2, H_f.cols()).noalias() = Hf_point;
				H_x.block(row, map_hx[clone_Ii], 2, clone_Ii->getSize()).noalias() = Hx_point;

				for (size_t i = 0; i < dpfg_dx_order.size(); i++)
				{
					H_x.block(row, map_hx[dpfg_dx_order.at(i)], 2, dpfg_dx_order.at(i)->getSize()).noalias() += Hx_point_extra.at(i);
				}

				if (state_ptr->options.do_calib_camera_pose)
				{
					H_x.block(row, map_hx[calibration], 2, calibration->getSize()).noalias() += dz_dpfc * dpfc_dcalib;
				}

				if (state_ptr->options.do_calib_camera_intrinsics)
				{
					H_x.block(row, map_hx[distortion], 2, distortion->getSize()) = Hx_dist;
				}
			}

			if (state_ptr->options.use_huber)
			{
				Eigen::VectorXd res_block = res.block(row, 0, block_rows, 1);
				double loss = res_block.norm();
				double hw = loss < setting_huber_th ? 1 : setting_huber_th / loss;

				if (hw < 1) hw = sqrt(hw);

				res.block(row, 0, block_rows, 1) *= hw;
				H_x.block(row, 0, block_rows, H_x.cols()) *= hw;
				H_f.block(row, 0, block_rows, H_f.cols()) *= hw;
			}

			row += block_rows;
		}
	}
}

void updaterHelper::nullspaceProjectInplace(Eigen::MatrixXd &H_f, Eigen::MatrixXd &H_x, Eigen::VectorXd &res)
{
	Eigen::JacobiRotation<double> Ho_GR_temp;
  
	for (int n = 0; n < H_f.cols(); ++n)
	{
		for (int m = (int)H_f.rows() - 1; m > n; m--)
		{
			Ho_GR_temp.makeGivens(H_f(m - 1, n), H_f(m, n));
			(H_f.block(m - 1, n, 2, H_f.cols() - n)).applyOnTheLeft(0, 1, Ho_GR_temp.adjoint());
			(H_x.block(m - 1, 0, 2, H_x.cols())).applyOnTheLeft(0, 1, Ho_GR_temp.adjoint());
			(res.block(m - 1, 0, 2, 1)).applyOnTheLeft(0, 1, Ho_GR_temp.adjoint());
		}
	}

	H_x = H_x.block(H_f.cols(), 0, H_x.rows() - H_f.cols(), H_x.cols()).eval();
	res = res.block(H_f.cols(), 0, res.rows() - H_f.cols(), res.cols()).eval();

	assert(H_x.rows() == res.rows());
}

void updaterHelper::measurementCompressInplace(Eigen::MatrixXd &H_x, Eigen::VectorXd &res)
{
	if (H_x.rows() <= H_x.cols())
		return;

	Eigen::JacobiRotation<double> Ho_GR_temp;
	for (int n = 0; n < H_x.cols(); n++)
	{
		for (int m = (int)H_x.rows() - 1; m > n; m--)
		{
			Ho_GR_temp.makeGivens(H_x(m - 1, n), H_x(m, n));

			(H_x.block(m - 1, n, 2, H_x.cols() - n)).applyOnTheLeft(0, 1, Ho_GR_temp.adjoint());
			(res.block(m - 1, 0, 2, 1)).applyOnTheLeft(0, 1, Ho_GR_temp.adjoint());
		}
	}

	int r = std::min(H_x.rows(), H_x.cols());

	assert(r <= H_x.rows());
	H_x.conservativeResize(r, H_x.cols());
	res.conservativeResize(r, res.cols());
}
