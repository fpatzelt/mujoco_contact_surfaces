/**
 *
 * Parts of the code used to evaluate the contact surfaces
 * (namely the methods evaluateContactSurfaces and passive_cb)
 * are based on code of Drake which is licensed as follows:
 *
 * All components of Drake are licensed under the BSD 3-Clause License
 * shown below. Where noted in the source code, some portions may
 * be subject to other permissive, non-viral licenses.
 *
 * Copyright 2012-2022 Robot Locomotion Group @ CSAIL
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.  Redistributions
 * in binary form must reproduce the above copyright notice, this list of
 * conditions and the following disclaimer in the documentation and/or
 * other materials provided with the distribution.  Neither the name of
 * the Massachusetts Institute of Technology nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 * The rest is licensed under the following license:
 *
 * Software License Agreement (BSD 3-Clause License)
 *
 *  Copyright (c) 2022, Bielefeld University
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *
 *   * Neither the name of Bielefeld University nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

/* Authors: Florian Patzelt*/

#include <mujoco_contact_surfaces/mujoco_contact_surfaces_plugin.h>

#include <pluginlib/class_list_macros.h>

namespace mujoco_contact_surfaces {

namespace {
// Map from data pointers to instances of the plugin. This instance is used in callback wrappers that are called from
// mujoco in order to find the correct plugin instance to use.
std::map<const mjData *, MujocoContactSurfacesPlugin *> instance_map;
mjfCollision defaultCollisionFunctions[mjNGEOMTYPES][mjNGEOMTYPES];
mjfGeneric default_passive_cb;

int collision_cb_wrapper(const mjModel *m, const mjData *d, mjContact *con, int g1, int g2, mjtNum margin)
{
	return instance_map[d]->collision_cb(m, d, con, g1, g2, margin);
}

void passive_cb_wrapper(const mjModel *m, mjData *d)
{
	if (default_passive_cb != NULL) {
		default_passive_cb(m, d);
	}
	if (instance_map.find(d) != instance_map.end()) {
		instance_map[d]->passive_cb(m, d);
	}
}

SpatialVelocity<double> getGeomVelocity(int id, const mjModel *m, const mjData *d)
{
	mjtNum res[6];
	mj_objectVelocity(m, d, mjOBJ_GEOM, id, res, 0);
	Vector3<double> w = Vector3<double>(res[0], res[1], res[2]);
	Vector3<double> v = Vector3<double>(res[3], res[4], res[5]);
	return SpatialVelocity<double>(w, v);
}

RigidTransform<double> getGeomPose(int id, const mjData *d)
{
	Eigen::Matrix3d m_;
	m_ << d->geom_xmat[9 * id + 0], d->geom_xmat[9 * id + 1], d->geom_xmat[9 * id + 2], d->geom_xmat[9 * id + 3],
	    d->geom_xmat[9 * id + 4], d->geom_xmat[9 * id + 5], d->geom_xmat[9 * id + 6], d->geom_xmat[9 * id + 7],
	    d->geom_xmat[9 * id + 8];

	return RigidTransform<double>(
	    RotationMatrix<double>(m_),
	    Vector3<double>(d->geom_xpos[3 * id], d->geom_xpos[3 * id + 1], d->geom_xpos[3 * id + 2]));
}

double calcCombinedElasticModulus(double E_A, double E_B)
{
	const double kInf = std::numeric_limits<double>::infinity();
	if (E_A == kInf)
		return E_B;
	if (E_B == kInf)
		return E_A;
	return E_A * E_B / (E_A + E_B);
}

double calcCombinedDissipation(ContactProperties *cp1, ContactProperties *cp2)
{
	const double kInf = std::numeric_limits<double>::infinity();

	const double E_A   = cp1->hydroelastic_modulus;
	const double E_B   = cp2->hydroelastic_modulus;
	const double d_A   = cp1->dissipation;
	const double d_B   = cp2->dissipation;
	const double Estar = calcCombinedElasticModulus(E_A, E_B);

	// Both bodies are rigid. We simply return the arithmetic average.
	if (Estar == kInf)
		return 0.5 * (d_A + d_B);

	// At least one body is soft.
	double d_star = 0;
	if (E_A != kInf)
		d_star += Estar / E_A * d_A;
	if (E_B != kInf)
		d_star += Estar / E_B * d_B;
	return d_star;
}

Vector3<double> CalcCentroidOfEnclosedVolume(const TriangleSurfaceMesh<double> &surface_mesh)
{
	double six_total_volume = 0;
	Vector3<double> centroid(0, 0, 0);

	// For convenience we tetrahedralize each triangle (p,q,r) about the
	// origin o = (0,0,0) in the mesh's frame. The centroid is then the
	// signed volume weighted sum of the centroids of the tetrahedra. The
	// choice of o is arbitrary and need not be on the interior of surface_mesh.
	// For efficiency we compute 6*Vi for the signed volume Vi of the ith element
	// and 6*V for the signed volume V of the region.
	for (const auto &tri : surface_mesh.triangles()) {
		const Vector3<double> &p = surface_mesh.vertex(tri.vertex(0));
		const Vector3<double> &q = surface_mesh.vertex(tri.vertex(1));
		const Vector3<double> &r = surface_mesh.vertex(tri.vertex(2));

		// 6 * signed volume of (p,q,r,o)
		const double six_volume = (p).cross(q).dot(r);
		six_total_volume += six_volume;

		// Centroid of tetrahedron (p,q,r,o) is (p + q + r + o) / 4.
		// We factor out the division for added precision and accuracy.
		centroid += six_volume * (p + q + r);
	}

	return centroid / (4 * six_total_volume);
}

} // namespace

MujocoContactSurfacesPlugin::~MujocoContactSurfacesPlugin()
{
	ROS_INFO_STREAM_NAMED("mujoco_contact_surfaces", "Shutting down mujoco_contact_surfaces plugin ...");
	geomCollisions.clear();
	contactProperties.clear();
	instance_map.erase(d_.get());
	if (instance_map.empty()) {
		for (int i = 0; i < mjNGEOMTYPES; ++i) {
			for (int j = 0; j < mjNGEOMTYPES; ++j) {
				if (mjCOLLISIONFUNC[i][j] == collision_cb_wrapper) {
					mjCOLLISIONFUNC[i][j] = defaultCollisionFunctions[i][j];
				}
			}
		}
	}
	delete this;
}

bool MujocoContactSurfacesPlugin::load(mjModelPtr m, mjDataPtr d)
{
	ROS_INFO_STREAM_NAMED("mujoco_contact_surfaces", "Loading mujoco_contact_surfaces plugin ...");
	d_                    = d;
	m_                    = m;
	instance_map[d.get()] = this;
	if (instance_map.size() == 1) {
		initCollisionFunction();
	}
	parseMujocoCustomFields(m.get());

	// TODO: This should not be hardcoded. Add tactile sensor
	int id = mj_name2id(m.get(), mjOBJ_GEOM, "myrmex_foam");
	if (id >= 0) {
		TactileSensor *ts = new TactileSensor();
		ts->geomID        = id;
		ts->geomName      = "myrmex_foam";
		ts->resolution    = 0.025;
		// for (double x = -0.19; x < 0.2; x = x + 0.02) {
		// 	for (double y = -0.19; y < 0.2; y = y + 0.02) {
		// 		ts->cellLocations.push_back(Vector3<double>(x, y, 0));
		// 	}
		// }
		// ts->cellLocations = { Vector3<double>(0, 0, 0) };
		tactileSensors.push_back(ts);
	}

	ROS_INFO_NAMED("mujoco_contact_surfaces", "Loaded mujoco_contact_surfaces");
	return true;
}

void MujocoContactSurfacesPlugin::update() {}

void MujocoContactSurfacesPlugin::reset()
{
	geomCollisions.clear();
	// contactProperties.clear();
	// instance_map.erase(d_.get());
}

int MujocoContactSurfacesPlugin::collision_cb(const mjModel *m, const mjData *d, mjContact *con, int g1, int g2,
                                              mjtNum margin)
{
	int t1                 = m->geom_type[g1];
	int t2                 = m->geom_type[g2];
	ContactProperties *cp1 = contactProperties[g1];
	ContactProperties *cp2 = contactProperties[g2];
	if (cp1 == NULL or cp2 == NULL or (cp1->contact_type == RIGID and cp2->contact_type == RIGID)) {
		return defaultCollisionFunctions[t1][t2](m, d, con, g1, g2, margin);
	}

	std::unique_ptr<ContactSurface<double>> s;
	if (cp1->contact_type == SOFT and cp2->contact_type == SOFT) {
		RigidTransform<double> p1 = getGeomPose(g1, d);
		RigidTransform<double> p2 = getGeomPose(g2, d);
		s = ComputeContactSurfaceFromCompliantVolumes<double>(cp1->drake_id, *cp1->pf, *cp1->bvh_v, p1, cp2->drake_id,
		                                                      *cp2->pf, *cp2->bvh_v, p2,
		                                                      hydroelastic_contact_representation);
	} else {
		if (cp1->contact_type == RIGID) {
			std::swap(cp1, cp2);
			std::swap(g1, g2);
		}
		RigidTransform<double> p1 = getGeomPose(g1, d);
		RigidTransform<double> p2 = getGeomPose(g2, d);
		if (m->geom_type[g2] == mjGEOM_PLANE) {
			s = ComputeContactSurfaceFromSoftVolumeRigidHalfSpace(cp1->drake_id, *cp1->pf, *cp1->bvh_v, p1, cp2->drake_id,
			                                                      p2, hydroelastic_contact_representation);
		} else {
			s = ComputeContactSurfaceFromSoftVolumeRigidSurface<double>(cp1->drake_id, *cp1->pf, *cp1->bvh_v, p1,
			                                                            cp2->drake_id, *cp2->sm, *cp2->bvh_s, p2,
			                                                            hydroelastic_contact_representation);
		}
	}

	if (s != NULL) {
		if (cp2->drake_id < cp1->drake_id) {
			std::swap(cp1, cp2);
			std::swap(g1, g2);
		}
		GeomCollision *gc = new GeomCollision(g1, g2, *s.get());
		evaluateContactSurface(m, d, gc);
		geomCollisions.push_back(gc);
	}

	return 0;
}

void MujocoContactSurfacesPlugin::evaluateContactSurface(const mjModel *m, const mjData *d, GeomCollision *gc)
{
	ContactSurface<double> *s = &gc->s;
	int g1                    = gc->g1;
	int g2                    = gc->g2;
	ContactProperties *cp1    = contactProperties[g1];
	ContactProperties *cp2    = contactProperties[g2];
	const double dissipation  = calcCombinedDissipation(cp1, cp2);
	// const double stiction_tolerance = 1.0e-4;
	// const double relative_tolerance = 1.0e-2;
	for (int face = 0; face < s->num_faces(); ++face) {
		const double &Ae = s->area(face); // Face element area.

		// We found out that the hydroelastic query might report
		// infinitesimally small triangles (consider for instance an initial
		// condition that perfectly places an object at zero distance from the
		// ground.) While the area of zero sized triangles is not a problem by
		// itself, the badly computed normal on these triangles leads to
		// problems when computing the contact Jacobians (since we need to
		// obtain an orthonormal basis based on that normal.)
		// We therefore ignore infinitesimally small triangles. The tolerance
		// below is somehow arbitrary and could possibly be tightened.
		if (Ae > 1.0e-14) {
			// From ContactSurface's documentation: The normal of each face is
			// guaranteed to point "out of" N and "into" M.
			const Vector3<double> &nhat_W = s->face_normal(face);

			// One dimensional pressure gradient (in Pa/m). Unlike [Masterjohn
			// et al. 2021], for convenience we define both pressure gradients
			// to be positive in the direction "into" the bodies. Therefore,
			// we use the minus sign for gN.
			// [Masterjohn et al., 2021] Discrete Approximation of Pressure
			// Field Contact Patches.
			const double gM =
			    s->HasGradE_M() ? s->EvaluateGradE_M_W(face).dot(nhat_W) : double(std::numeric_limits<double>::infinity());
			const double gN = s->HasGradE_N() ? -s->EvaluateGradE_N_W(face).dot(nhat_W) :
			                                    double(std::numeric_limits<double>::infinity());

			constexpr double kGradientEpsilon = 1.0e-14;
			if (gM < kGradientEpsilon || gN < kGradientEpsilon) {
				// Mathematically g = gN*gM/(gN+gM) and therefore g = 0 when
				// either gradient on one of the bodies is zero. A zero gradient
				// means there is no contact constraint, and therefore we
				// ignore it to avoid numerical problems in the discrete solver.
				continue;
			}

			// Effective hydroelastic pressure gradient g result of
			// compliant-compliant interaction, see [Masterjohn et al., 2021].
			// The expression below is mathematically equivalent to g =
			// gN*gM/(gN+gM) but it has the advantage of also being valid if
			// one of the gradients is infinity.
			const double g = 1.0 / (1.0 / gM + 1.0 / gN);

			// Position of quadrature point Q in the world frame (since mesh_W
			// is measured and expressed in W).
			const Vector3<double> &p_WQ = s->centroid(face);
			// For a triangle, its centroid has the fixed barycentric
			// coordinates independent of the shape of the triangle. Using
			// barycentric coordinates to evaluate field value could be
			// faster than using Cartesian coordiantes, especially if the
			// TriangleSurfaceMeshFieldLinear<> does not store gradients and
			// has to solve linear equations to convert Cartesian to
			// barycentric coordinates.
			const Vector3<double> tri_centroid_barycentric(1 / 3., 1 / 3., 1 / 3.);
			// Pressure at the quadrature point.
			const double p0 = s->is_triangle() ? s->tri_e_MN().Evaluate(face, tri_centroid_barycentric) :
			                                     s->poly_e_MN().EvaluateCartesian(face, p_WQ);

			// Force contribution by this quadrature point.
			const double fn0 = Ae * p0;

			// Effective compliance in the normal direction for the given
			// discrete patch, refer to [Masterjohn et al., 2021] for details.
			// [Masterjohn, 2021] Masterjohn J., Guoy D., Shepherd J. and Castro
			// A., 2021. Discrete Approximation of Pressure Field Contact Patches.
			// Available at https://arxiv.org/abs/2110.04157.
			const double k = Ae * g;

			gc->pointCollisions.push_back({ p_WQ, nhat_W, fn0, k, dissipation, face });
		}
	}
}

void MujocoContactSurfacesPlugin::passive_cb(const mjModel *m, mjData *d)
{
	if (visualizeContactSurfaces) {
		// reset the visualized geoms
		n_vGeom = 0;
		running_scale = 0.9 * running_scale + 0.1 * current_scale;
		current_scale = 0.;
	}

	const double stiction_tolerance = 1.0e-4; // TODO should this be hardcoded?
	const double relative_tolerance = 1.0e-2;
	for (GeomCollision *gc : geomCollisions) {
		int g1                 = gc->g1;
		int g2                 = gc->g2;
		ContactProperties *cp1 = contactProperties[g1];
		ContactProperties *cp2 = contactProperties[g2];
		const CoulombFriction<double> &geometryM_friction =
		    CoulombFriction<double>{ cp1->static_friction, cp1->dynamic_friction };
		const CoulombFriction<double> &geometryN_friction =
		    CoulombFriction<double>{ cp2->static_friction, cp2->dynamic_friction };
		const CoulombFriction<double> combined_friction =
		    CalcContactFrictionFromSurfaceProperties(geometryM_friction, geometryN_friction);
		const double mu_coulomb = combined_friction.dynamic_friction();
		for (PointCollision pc : gc->pointCollisions) {
			const RigidTransform<double> &X_WA  = getGeomPose(g1, d);
			const RigidTransform<double> &X_WB  = getGeomPose(g2, d);
			const SpatialVelocity<double> &V_WA = getGeomVelocity(g1, m, d);
			const SpatialVelocity<double> &V_WB = getGeomVelocity(g2, m, d);

			const Vector3<double> p_AoAq_W      = pc.p - X_WA.translation();
			const SpatialVelocity<double> V_WAq = V_WA.Shift(p_AoAq_W);

			// Next compute the spatial velocity of Body B at Bq.
			const Vector3<double> p_BoBq_W      = pc.p - X_WB.translation();
			const SpatialVelocity<double> V_WBq = V_WB.Shift(p_BoBq_W);

			// Finally compute the relative velocity of Frame Aq relative to Frame Bq,
			// expressed in the world frame, and then the translational component of this
			// velocity.
			const SpatialVelocity<double> V_BqAq_W = V_WAq - V_WBq;
			const Vector3<double> &v_BqAq_W        = V_BqAq_W.translational();

			// Get the velocity along the normal to the contact surface. Note that a
			// positive value indicates that bodies are separating at Q while a negative
			// value indicates that bodies are approaching at Q.
			const double vn_BqAq_W = v_BqAq_W.dot(pc.n);

			const double fn = std::max(0., 1. - pc.damping * vn_BqAq_W) * (pc.fn0 - 0.001 * pc.stiffness * vn_BqAq_W);

			const Vector3<double> vt   = v_BqAq_W - pc.n * vn_BqAq_W;
			double epsilon             = stiction_tolerance * relative_tolerance;
			epsilon                    = epsilon * epsilon;
			const double v_slip        = std::sqrt(vt.squaredNorm() + epsilon);
			const Vector3<double> that = vt / v_slip;
			double mu_regularized      = mu_coulomb;
			const double s             = v_slip / stiction_tolerance;
			if (s < 1) {
				mu_regularized = mu_coulomb * s * (2.0 - s);
			}

			const Vector3<double> f_slip = -mu_regularized * that * fn;

			const Vector3<double> f = f_slip + fn * pc.n;

			const mjtNum point[3]  = { pc.p[0], pc.p[1], pc.p[2] };
			const mjtNum torque[3] = { 0, 0, 0 };
			const mjtNum forceA[3] = { f[0], f[1], f[2] };
			const mjtNum forceB[3] = { -f[0], -f[1], -f[2] };
			mj_applyFT(m, d, forceA, torque, point, m->geom_bodyid[g1], d->qfrc_passive);
			mj_applyFT(m, d, forceB, torque, point, m->geom_bodyid[g2], d->qfrc_passive);

			// visualize collision force:
			// int id = contactProperties[g1]->contact_type == SOFT ? g1 : g2;
			// // ContactProperties *cp =
			// //     contactProperties[g1]->contact_type == SOFT ? contactProperties[g1] : contactProperties[g2];
			// if (m->geom_type[id] == mjGEOM_BOX) {
			// 	// ROS_INFO_STREAM_NAMED("mujoco_contact_surfaces", "fn: " << fn);
			// 	const float rgba[4] = { std::min(1. - fn, 1.), 0, std::max(fn, 0.0), 0.8 };
			// 	mjtNum pos[3];
			// 	mjtNum size[3];
			// 	for (int i = 0; i < 3; ++i) {
			// 		pos[i]  = d->geom_xpos[3 * id + i];
			// 		size[i] = m->geom_size[3 * id + i];
			// 	}
			// 	mjtNum rot[9];
			// 	for (int i = 0; i < 9; ++i) {
			// 		rot[i] = d->geom_xmat[9 * id + i];
			// 	}
			// 	if (n_vGeom == MAX_VGEOM) {
			// 		break;
			// 	}
			// 	mjvGeom *g = vGeoms + n_vGeom++;
			// 	mjv_initGeom(g, mjGEOM_BOX, size, pos, rot, rgba);
			// }
			// visualize
			if (visualizeContactSurfaces) {
				current_scale = std::max(current_scale, std::abs(fn));
				if (gc->s.is_triangle()) {
					visualizeMeshElement(pc.face, gc->s.tri_mesh_W(), fn);
				} else {
					visualizeMeshElement(pc.face, gc->s.poly_mesh_W(), fn);
				}
			}
		}
		// Tactile sensors
		// First approach: Interate over sensor cells
		// for (TactileSensor *ts : tactileSensors) {
		// 	int id = ts->geomID;
		// 	if (g1 == ts->geomID or g2 == ts->geomID) {
		// 		ContactSurface<double> s = gc->s;
		// 		auto mesh                = s.tri_mesh_W();
		// 		for (Vector3<double> p : ts->cellLocations) {
		// 			for (int t = 0; t < mesh.num_elements(); ++t) {
		// 				Vector3<double> b = mesh.CalcBarycentric(p, t).normalized();

		// 				if (b.minCoeff() >= 0) {
		// 					ROS_INFO_STREAM_NAMED("mujoco_contact_surfaces", "Bary: " << b);
		// 					const float rgba[4] = { 1, 0, 0, 0.8 };
		// 					mjtNum pos[3];
		// 					mjtNum size[3] = { 0.01, 0.01, 0.003 };
		// 					for (int i = 0; i < 3; ++i) {
		// 						pos[i] = d->geom_xpos[3 * id + i] + p[i];
		// 					}
		// 					mjtNum rot[9];
		// 					for (int i = 0; i < 9; ++i) {
		// 						rot[i] = d->geom_xmat[9 * id + i];
		// 					}
		// 					if (n_vGeom == MAX_VGEOM) {
		// 						break;
		// 					}
		// 					mjvGeom *g = vGeoms + n_vGeom++;
		// 					mjv_initGeom(g, mjGEOM_BOX, size, pos, rot, rgba);
		// 					break;
		// 				}
		// 			}
		// 		}
		// 	}

		// second approach: iterate over vertices
		// for (TactileSensor *ts : tactileSensors) {
		// 	int id        = ts->geomID;
		// 	double height = m->geom_size[3 * id + 2];
		// 	double h2     = height * height;
		// 	if (g1 == ts->geomID or g2 == ts->geomID) {
		// 		ContactSurface<double> s = gc->s;
		// 		auto mesh                = s.tri_mesh_W();
		// 		for (int t = 0; t < mesh.num_elements(); ++t) {
		// 			auto tri                  = mesh.element(t);
		// 			const Vector3<double> &v0 = mesh.vertex(tri.vertex(0));
		// 			const Vector3<double> &v1 = mesh.vertex(tri.vertex(1));
		// 			const Vector3<double> &v2 = mesh.vertex(tri.vertex(2));
		// 			const Vector3<double> &c  = mesh.element_centroid(t);

		// 			double d0 = (v0 - c).norm();
		// 			double d1 = (v1 - c).norm();
		// 			double d2 = (v2 - c).norm();

		// 			double dist = std::max({ d0, d1, d2 });
		// 			dist        = std::sqrt(dist * dist + h2);
		// 			Eigen::Matrix<double, 3, 2> A;
		// 			A.col(0) << v1 - v0;
		// 			A.col(1) << v2 - v0;
		// 			auto h = A.colPivHouseholderQr();

		// 			for (Vector3<double> p0 : ts->cellLocations) {
		// 				Vector3<double> p;
		// 				for (int i = 0; i < 3; ++i) {
		// 					p[i] = d->geom_xpos[3 * id + i] + p0[i];
		// 				}
		// 				if ((p - c).norm() <= dist) {
		// 					Vector2<double> solution = h.solve(p - v0);
		// 					const double &b1         = solution(0);
		// 					const double &b2         = solution(1);
		// 					const double b0          = 1. - b1 - b2;
		// 					const Vector3<double> b  = Vector3<double>(b0, b1, b2).normalized();
		// 					if (b.minCoeff() >= 0) {
		// 						// ROS_INFO_STREAM_NAMED("mujoco_contact_surfaces", "Bary: " << b);
		// 						const float rgba[4] = { 1, 0, 0, 0.8 };
		// 						mjtNum pos[3]       = { p[0], p[1], p[2] };
		// 						mjtNum size[3]      = { 0.01, 0.01, 0.003 };
		// 						mjtNum rot[9];
		// 						for (int i = 0; i < 9; ++i) {
		// 							rot[i] = d->geom_xmat[9 * id + i];
		// 						}
		// 						if (n_vGeom < MAX_VGEOM) {
		// 							mjvGeom *g = vGeoms + n_vGeom++;
		// 							mjv_initGeom(g, mjGEOM_BOX, size, pos, rot, rgba);
		// 						}
		// 					}
		// 				}
		// 			}
		// 		}
		// 	}
		// }

		// third approach: caching
		// for (TactileSensor *ts : tactileSensors) {
		// 	int id        = ts->geomID;
		// 	double height = m->geom_size[3 * id + 2];
		// 	double h2     = height * height;
		// 	if (g1 == ts->geomID or g2 == ts->geomID) {
		// 		ContactSurface<double> s = gc->s;
		// 		auto mesh                = s.tri_mesh_W();

		// 		// prepare caches
		// 		const int m = ts->cellLocations.size();
		// 		const int n = mesh.num_elements();
		// 		Eigen::MatrixX3d centroids(n, 3);
		// 		Eigen::MatrixX3d cells(m, 3);
		// 		Eigen::VectorXd distance_thresholds(n);
		// 		std::map<int, Eigen::ColPivHouseholderQR<Eigen::Matrix<double, 3, 2>>> cph_cache = {};
		// 		std::vector<Vector3<double>> v0s;
		// 		std::vector<Vector3<double>> v1s;
		// 		std::vector<Vector3<double>> v2s;

		// 		for (int i = 0; i < m; ++i) {
		// 			Vector3<double> p0 = ts->cellLocations[i];
		// 			Vector3<double> p;
		// 			for (int j = 0; j < 3; ++j) {
		// 				p[j] = d->geom_xpos[3 * id + j] + p0[j];
		// 			}
		// 			cells.row(i) << p.transpose();
		// 		}

		// 		auto time0 = high_resolution_clock::now();

		// 		for (int t = 0; t < n; ++t) {
		// 			auto tri                  = mesh.element(t);
		// 			const Vector3<double> &v0 = mesh.vertex(tri.vertex(0));
		// 			const Vector3<double> &v1 = mesh.vertex(tri.vertex(1));
		// 			const Vector3<double> &v2 = mesh.vertex(tri.vertex(2));
		// 			const Vector3<double> &c  = mesh.element_centroid(t);
		// 			v0s.push_back(v0);
		// 			v1s.push_back(v1);
		// 			v2s.push_back(v2);
		// 			centroids.row(t) << c.transpose();

		// 			double d0 = (v0 - c).norm();
		// 			double d1 = (v1 - c).norm();
		// 			double d2 = (v2 - c).norm();

		// 			double dist            = std::max({ d0, d1, d2 });
		// 			dist                   = std::sqrt(dist * dist + h2);
		// 			distance_thresholds[t] = dist;
		// 		}
		// 		auto time1 = high_resolution_clock::now();
		// 		ROS_INFO_STREAM_NAMED("mujoco_contact_surfaces", "Caching time: " << duration_cast<microseconds>(time1 -
		// time0).count()); 		Eigen::MatrixXd distances(n,m);
		// 		//auto distances = centroids * cells.transpose();//
		// 		distances = (((-2 * centroids * cells.transpose()).colwise() +
		// centroids.transpose().colwise().squaredNorm().transpose()).rowwise() +
		// cells.transpose().colwise().squaredNorm()).cwiseSqrt().colwise() - distance_thresholds;
		// 		//distances = distances.colwise() - distance_thresholds;
		// 		auto time1_ = high_resolution_clock::now();
		// 		ROS_INFO_STREAM_NAMED("mujoco_contact_surfaces", "Distances: " << duration_cast<microseconds>(time1_ -
		// time1).count());

		// 		for (int l = 0; l < m; ++l) {
		// 			Vector3<double> p = cells.row(l);
		// 			auto time10 = high_resolution_clock::now();
		// 			const Eigen::VectorXd dists = distances.col(l);
		// 			auto time11 = high_resolution_clock::now();
		// 			for (int t = 0; t < n; ++t) {
		// 				if (dists[t] <= 0) {
		// 					if (cph_cache.count(t) == 0) {
		// 						Eigen::Matrix<double, 3, 2> A;
		// 						A.col(0) << v1s[t] - v0s[t];
		// 						A.col(1) << v2s[t] - v0s[t];
		// 						cph_cache[t] = A.colPivHouseholderQr();
		// 					}

		// 					Vector2<double> solution = cph_cache[t].solve(p - v0s[t]);
		// 					const double &b1         = solution(0);
		// 					const double &b2         = solution(1);
		// 					const double b0          = 1. - b1 - b2;
		// 					const Vector3<double> b  = Vector3<double>(b0, b1, b2).normalized();
		// 					if (b.minCoeff() >= 0) {
		// 						// ROS_INFO_STREAM_NAMED("mujoco_contact_surfaces", "Bary: " << b.transpose());
		// 						const float rgba[4] = { 1, 0, 0, 0.8 };
		// 						mjtNum pos[3]       = { p[0], p[1], p[2] };
		// 						mjtNum size[3]      = { 0.01, 0.01, 0.003 };
		// 						mjtNum rot[9];
		// 						for (int i = 0; i < 9; ++i) {
		// 							rot[i] = d->geom_xmat[9 * id + i];
		// 						}
		// 						if (n_vGeom < MAX_VGEOM) {
		// 							mjvGeom *g = vGeoms + n_vGeom++;
		// 							mjv_initGeom(g, mjGEOM_BOX, size, pos, rot, rgba);
		// 						}
		// 					}
		// 				}
		// 			}
		// 			auto time12 = high_resolution_clock::now();
		// 			//ROS_INFO_STREAM_NAMED("mujoco_contact_surfaces", "Loop1: " << duration_cast<microseconds>(time11 -
		// time10).count() << " " << duration_cast<microseconds>(time12 - time11).count());
		// 		}
		// 		auto time2 = high_resolution_clock::now();
		// 		ROS_INFO_STREAM_NAMED("mujoco_contact_surfaces", "TS time: " << duration_cast<microseconds>(time2 -
		// time1_).count());
		// 	}
		// }

		for (TactileSensor *ts : tactileSensors) {
			int id        = ts->geomID;
			double height = m->geom_size[3 * id + 2];
			double h2     = height * height;
			if (g1 == ts->geomID or g2 == ts->geomID) {
				ContactSurface<double> s = gc->s;
				auto mesh                = s.tri_mesh_W();
				double xs                = m->geom_size[3 * id];
				double ys                = m->geom_size[3 * id + 1];
				double zs                = m->geom_size[3 * id + 2];
				double res               = ts->resolution;
				int cx                   = (int)(2 * xs / res) + 1;
				int cy                   = (int)(2 * ys / res) + 1;

				std::vector<Eigen::Vector4d> cpoints[cx][cy];

				mjtNum size[3] = { res / 2, res / 2, 0.001 };
				mjtNum rot[9];
				for (int i = 0; i < 9; ++i) {
					rot[i] = d->geom_xmat[9 * id + i];
				}

				// prepare caches
				const int n = mesh.num_elements();
				const int m = mesh.num_vertices();

				// bool used_vertices[m];
				// std::fill_n(used_vertices, m, false);

				// get geom transformation
				Eigen::Matrix4d M, Minv, Mback;
				M << d->geom_xmat[9 * id + 0], d->geom_xmat[9 * id + 1], d->geom_xmat[9 * id + 2], d->geom_xpos[3 * id],
				    d->geom_xmat[9 * id + 3], d->geom_xmat[9 * id + 4], d->geom_xmat[9 * id + 5], d->geom_xpos[3 * id + 1],
				    d->geom_xmat[9 * id + 6], d->geom_xmat[9 * id + 7], d->geom_xmat[9 * id + 8], d->geom_xpos[3 * id + 2],
				    0, 0, 0, 1;

				Minv = M.inverse();
				Minv.col(3) << Minv.col(3) + Eigen::Vector4d(xs, ys, -zs, 0);

				Mback = Minv.inverse();
				int tris[cx][cy];
				Vector3<double> barys[cx][cy];
				double dists[cx][cy];
				for (int i = 0; i < cx; ++i) {
					for (int j = 0; j < cy; ++j) {
						dists[i][j] = res;
					}
				}

				for (int t = 0; t < n; ++t) {
					// tpoints.push_back(mesh.element_centroid(t));
					// project points onto 2d sensor plane
					std::vector<Eigen::Vector2d> tpoints = {};
					auto element                         = mesh.element(t);
					for (int i = 0; i < element.num_vertices(); ++i) {
						int v                     = element.vertex(i);
						const Vector3<double> &vp = mesh.vertex(v);
						const Eigen::Vector4d vpe(vp[0], vp[1], vp[2], 1);
						Eigen::Vector4d vpp = Minv * vpe;

						// if (!used_vertices[v]) {
						// 	used_vertices[v] = true;
						tpoints.push_back(Eigen::Vector2d(vpp[0], vpp[1]));
						// }
					}
					int s0 = (int)((tpoints[1] - tpoints[0]).norm() / res * 2.) + 1;
					int s1 = (int)((tpoints[2] - tpoints[0]).norm() / res * 2.) + 1;
					int s2 = (int)((tpoints[2] - tpoints[0]).norm() / res * 2.) + 1;
					int s  = std::max(s0, s1);
					// ROS_INFO_STREAM_NAMED("mujoco_contact_surfaces", s);
					// s = s2 = 3;
					for (double a = 0; a <= 1; a += 1. / s) {
						for (double b = 0; b <= 1; b += 1. / s2) {
							const Vector3<double> bary(a, (1 - a) * (1 - b), (1 - a) * b);
							// bary.normalize();
							Eigen::Vector2d p = bary[0] * tpoints[0] + bary[1] * tpoints[1] + bary[2] * tpoints[2];
							if (p[0] > 0 && p[0] < 2 * xs && p[1] > 0 && p[1] < 2 * ys) {
								int x = (int)std::floor(p[0] / res);
								int y = (int)std::floor(p[1] / res);
								Eigen::Vector2d c(x * res + res / 2., y * res + res / 2.);
								double dist = (c - p).norm();
								if (dist < dists[x][y]) {
									dists[x][y] = dist;
									barys[x][y] = bary;
									tris[x][y]  = t;
								}
								// Eigen::Vector4d dp = Mback * Eigen::Vector4d(x * res + res / 2, y * res + res / 2, 0, 1);
								// mjtNum pos[3]      = { dp[0], dp[1], dp[2] };
								// if (n_vGeom < MAX_VGEOM) {
								// 	mjvGeom *g = vGeoms + n_vGeom++;
								// 	mjv_initGeom(g, mjGEOM_BOX, size, pos, rot, rgba);
								// }
							}
						}
					}
				}

				for (int x = 0; x < cx; ++x) {
					for (int y = 0; y < cy; ++y) {
						if (dists[x][y] < res) {
							int t               = tris[x][y];
							double p0           = s.tri_e_MN().Evaluate(t, barys[x][y]) * s.area(t);
							//ROS_INFO_STREAM_NAMED("mujoco_contact_surfaces", "P*A: " << p0);
							float ps            = std::min(1., (float)p0 / 0.025);
							const float rgba[4] = { ps, 0, (1.f - ps), 0.8 };
							Eigen::Vector4d dp  = Mback * Eigen::Vector4d(x * res + res / 2, y * res + res / 2, 0, 1);
							mjtNum pos[3]       = { dp[0], dp[1], dp[2] };
							if (n_vGeom < MAX_VGEOM) {
								mjvGeom *g = vGeoms + n_vGeom++;
								mjv_initGeom(g, mjGEOM_BOX, size, pos, rot, rgba);
							}
						}
					}
				}
			}
		}
	}

	geomCollisions.clear();
}



template <class T>
void MujocoContactSurfacesPlugin::visualizeMeshElement(int face, T mesh, double fn)
{
	if (n_vGeom < MAX_VGEOM) {
		const mjtNum size[3] = { 0.00015, 0.00015, 0.00015 };
		// int face                    = pc.face;
		float scale                 = std::min(std::abs(fn), running_scale) / running_scale;
		const float rgba[4]         = { scale, 0., 1.0f - scale, 0.8 };
		const Vector3<double> &p_WQ = mesh.element_centroid(face);
		const mjtNum pos[3]         = { p_WQ[0], p_WQ[1], p_WQ[2] };

		mjvGeom *g = vGeoms + n_vGeom++;
		mjv_initGeom(g, mjGEOM_SPHERE, size, pos, NULL, rgba);
		auto p              = mesh.element(face);
		Vector3<double> vp0 = mesh.vertex(p.vertex(p.num_vertices() - 1));
		Vector3<double> n   = mesh.face_normal(face);
		for (int v = 0; v < p.num_vertices(); ++v) {
			Vector3<double> vp1 = mesh.vertex(p.vertex(v));
			if (n_vGeom == MAX_VGEOM) {
				break;
			}
			mjvGeom *g = vGeoms + n_vGeom++;
			mjv_initGeom(g, mjGEOM_CYLINDER, NULL, NULL, NULL, rgba);
			mjv_makeConnector(g, mjGEOM_CYLINDER, 0.0001, vp0[0], vp0[1], vp0[2], vp1[0], vp1[1], vp1[2]);
			vp0 = vp1;
		}
	}
}

void MujocoContactSurfacesPlugin::passiveCallback(mjModelPtr model, mjDataPtr data)
{
	passive_cb(model.get(), data.get());
}

void MujocoContactSurfacesPlugin::initCollisionFunction()
{
	// #define SP std::placeholders
	// 	auto collision_function = std::bind(&MujocoContactSurfacesPlugin::collision_cb, this, SP::_1, SP::_2,
	// SP::_3, SP::_4, SP::_5, SP::_6);
	for (int i = 0; i < mjNGEOMTYPES; ++i) {
		for (int j = 0; j < mjNGEOMTYPES; ++j) {
			defaultCollisionFunctions[i][j] = mjCOLLISIONFUNC[i][j];
			// registerCollisionFunc(i, j, collision_function);
			registerCollisionFunc(i, j, collision_cb_wrapper);
		}
	}
}

void MujocoContactSurfacesPlugin::parseMujocoCustomFields(mjModel *m)
{
	// parse HydroelasticContactRepresentation
	int hcp_id = mj_name2id(m, mjOBJ_TEXT, (PREFIX + "HydroelasticContactRepresentation").c_str());
	if (hcp_id >= 0) {
		int hcp_adr = m->text_adr[hcp_id];
		if (hcp_adr >= 0) {
			int hcp_size = m->text_size[hcp_id];
			std::string hcp(&m->text_data[hcp_adr], hcp_size);
			if (hcp.find("kTriangle") != std::string::npos) {
				hydroelastic_contact_representation = HydroelasticContactRepresentation::kTriangle;
				ROS_INFO_STREAM_NAMED("mujoco_contact_surfaces", "Found HydroelasticContactRepresentation: kTriangle");
			} else if (hcp.find("kPolygon") != std::string::npos) {
				hydroelastic_contact_representation = HydroelasticContactRepresentation::kPolygon;
				ROS_INFO_STREAM_NAMED("mujoco_contact_surfaces", "Found HydroelasticContactRepresentation: kPolygon");
			} else {
				ROS_INFO_STREAM_NAMED("mujoco_contact_surfaces",
				                      "No HydroelasticContactRepresentation found. Using default: kPolygon");
			}
		}
	}
	// parse VisualizeSurfaces
	int vs_id = mj_name2id(m, mjOBJ_NUMERIC, (PREFIX + "VisualizeSurfaces").c_str());
	if (vs_id >= 0) {
		int vs_adr  = m->numeric_adr[vs_id];
		int vs_size = m->numeric_size[vs_adr];
		if (vs_size == 1) {
			visualizeContactSurfaces = m->numeric_data[vs_adr];
		}
		ROS_INFO_STREAM_NAMED("mujoco_contact_surfaces", "VisualizeContactSurfaces: " << visualizeContactSurfaces);
	}
	// parse geom contact properties
	for (int i = 0; i < m->nnumeric; ++i) {
		std::string full_name = mj_id2name(m, mjOBJ_NUMERIC, i);
		if (full_name.rfind(PREFIX, 0) == 0) {
			std::string s = full_name.substr(PREFIX.length());
			int id        = mj_name2id(m, mjOBJ_GEOM, s.c_str());
			if (id >= 0) {
				int adr  = m->numeric_adr[i];
				int size = m->numeric_size[i];
				if (adr >= 0 and size == 5) {
					double hydroElasticModulus = m->numeric_data[adr];
					double dissipation         = m->numeric_data[adr + 1];
					double resolutionHint      = m->numeric_data[adr + 2];
					double staticFriction      = m->numeric_data[adr + 3];
					double dynamicFriction     = m->numeric_data[adr + 4];
					ROS_INFO_STREAM_NAMED("mujoco_contact_surfaces",
					                      "Found geom '" << s << "' with properties: hM " << hydroElasticModulus << " dis "
					                                     << dissipation << " rH " << resolutionHint);
					ContactProperties *cp;
					int geomType = m->geom_type[id];

					switch (geomType) {
						case mjGEOM_PLANE:
							if (hydroElasticModulus > 0) {
								ROS_INFO_STREAM_NAMED("mujoco_contact_surfaces", "soft plane collision not implemented yet");
							} else {
								cp                    = new ContactProperties(id, s, RIGID, staticFriction, dynamicFriction);
								contactProperties[id] = cp;
							}
							break;
						case mjGEOM_HFIELD: // height field
							ROS_INFO_STREAM_NAMED("mujoco_contact_surfaces", "hfield collision not implemented yet");
							break;
						case mjGEOM_SPHERE: // sphere
						{
							Sphere *sphere = new Sphere(m->geom_size[3 * id]);
							if (hydroElasticModulus > 0) {
								VolumeMesh<double> *vm                    = new VolumeMesh<double>(MakeSphereVolumeMesh<double>(
                            *sphere, resolutionHint, TessellationStrategy::kSingleInteriorVertex));
								VolumeMeshFieldLinear<double, double> *pf = new VolumeMeshFieldLinear<double, double>(
								    MakeSpherePressureField<double>(*sphere, vm, hydroElasticModulus));
								Bvh<Obb, VolumeMesh<double>> *bvh = new Bvh<Obb, VolumeMesh<double>>(*vm);
								cp = new ContactProperties(id, s, SOFT, sphere, vm, pf, bvh, hydroElasticModulus, dissipation,
								                           staticFriction, dynamicFriction);
							} else {
								TriangleSurfaceMesh<double> *sm =
								    new TriangleSurfaceMesh<double>(MakeSphereSurfaceMesh<double>(*sphere, resolutionHint));
								Bvh<Obb, TriangleSurfaceMesh<double>> *bvh = new Bvh<Obb, TriangleSurfaceMesh<double>>(*sm);
								cp = new ContactProperties(id, s, RIGID, sphere, sm, bvh, staticFriction, dynamicFriction);
							}

							contactProperties[id] = cp;
							break;
						}
						case mjGEOM_CAPSULE: // capsule
							ROS_INFO_STREAM_NAMED("mujoco_contact_surfaces", "capsule collision not implemented yet");
							break;
						case mjGEOM_ELLIPSOID: // ellipsoid
							ROS_INFO_STREAM_NAMED("mujoco_contact_surfaces", "ellipsoid collision not implemented yet");
							break;
						case mjGEOM_CYLINDER: // cylinder
						{
							Cylinder *cylinder = new Cylinder(m->geom_size[3 * id], 2 * m->geom_size[3 * id + 1]);
							if (hydroElasticModulus > 0) {
								VolumeMesh<double> *vm;
								if (resolutionHint > 0) {
									vm = new VolumeMesh<double>(MakeCylinderVolumeMeshWithMa<double>(*cylinder, resolutionHint));
								}
								VolumeMeshFieldLinear<double, double> *pf = new VolumeMeshFieldLinear<double, double>(
								    MakeCylinderPressureField<double>(*cylinder, vm, hydroElasticModulus));
								Bvh<Obb, VolumeMesh<double>> *bvh = new Bvh<Obb, VolumeMesh<double>>(*vm);
								cp = new ContactProperties(id, s, SOFT, cylinder, vm, pf, bvh, hydroElasticModulus, dissipation,
								                           staticFriction, dynamicFriction);
							} else {
								TriangleSurfaceMesh<double> *sm =
								    new TriangleSurfaceMesh<double>(MakeCylinderSurfaceMesh<double>(*cylinder, resolutionHint));
								Bvh<Obb, TriangleSurfaceMesh<double>> *bvh = new Bvh<Obb, TriangleSurfaceMesh<double>>(*sm);
								cp = new ContactProperties(id, s, RIGID, cylinder, sm, bvh, staticFriction, dynamicFriction);
							}
							contactProperties[id] = cp;
							break;
						}
						case mjGEOM_BOX: // box
						{
							Box *box =
							    new Box(2 * m->geom_size[3 * id], 2 * m->geom_size[3 * id + 1], 2 * m->geom_size[3 * id + 2]);
							if (hydroElasticModulus > 0) {
								VolumeMesh<double> *vm;
								if (resolutionHint > 0) {
									vm = new VolumeMesh<double>(MakeBoxVolumeMesh<double>(*box, resolutionHint));
								} else {
									vm = new VolumeMesh<double>(MakeBoxVolumeMeshWithMa<double>(*box));
								}
								VolumeMeshFieldLinear<double, double> *pf = new VolumeMeshFieldLinear<double, double>(
								    MakeBoxPressureField<double>(*box, vm, hydroElasticModulus));
								Bvh<Obb, VolumeMesh<double>> *bvh = new Bvh<Obb, VolumeMesh<double>>(*vm);
								cp = new ContactProperties(id, s, SOFT, box, vm, pf, bvh, hydroElasticModulus, dissipation,
								                           staticFriction, dynamicFriction);
							} else {
								TriangleSurfaceMesh<double> *sm =
								    new TriangleSurfaceMesh<double>(MakeBoxSurfaceMesh<double>(*box, resolutionHint));
								Bvh<Obb, TriangleSurfaceMesh<double>> *bvh = new Bvh<Obb, TriangleSurfaceMesh<double>>(*sm);
								cp = new ContactProperties(id, s, RIGID, box, sm, bvh, staticFriction, dynamicFriction);
							}
							contactProperties[id] = cp;
							break;
						}
						case mjGEOM_MESH: // mesh
						{
							int geom_dataid = m->geom_dataid[id];
							if (geom_dataid >= 0) {
								int nv    = m->mesh_vertnum[geom_dataid];
								int nf    = m->mesh_facenum[geom_dataid];
								int v_adr = m->mesh_vertadr[geom_dataid];
								int f_adr = m->mesh_faceadr[geom_dataid];
								if (nv > 0 && nf > 0 && v_adr >= 0 && f_adr >= 0) {
									std::vector<SurfaceTriangle> triangles;
									std::vector<Vector3<double>> vertices;
									for (int i = 0; i < nv; ++i) {
										int v = v_adr + 3 * i;
										vertices.push_back(
										    Vector3<double>(m->mesh_vert[v], m->mesh_vert[v + 1], m->mesh_vert[v + 2]));
									}
									for (int i = 0; i < nf; ++i) {
										int f = f_adr + 3 * i;
										triangles.push_back(
										    SurfaceTriangle(m->mesh_face[f], m->mesh_face[f + 1], m->mesh_face[f + 2]));
									}
									TriangleSurfaceMesh<double> *sm =
									    new TriangleSurfaceMesh<double>(std::move(triangles), std::move(vertices));
									if (hydroElasticModulus > 0) {
										std::vector<Vector3<double>> volume_mesh_vertices(sm->vertices().begin(),
										                                                  sm->vertices().end());

										const Vector3<double> centroid = CalcCentroidOfEnclosedVolume(*sm);
										volume_mesh_vertices.push_back(centroid);

										const int centroid_index = volume_mesh_vertices.size() - 1;

										// The number of tetrahedra in the volume mesh is the same as the number of
										// triangles in the surface mesh.
										std::vector<VolumeElement> volume_mesh_elements;
										volume_mesh_elements.reserve(sm->num_elements());
										for (const SurfaceTriangle &e : sm->triangles()) {
											// Orient the tetrahedron such that it has positive signed volume (assuming
											// outward facing normal for the surface triangle).
											volume_mesh_elements.push_back(
											    { centroid_index, e.vertex(0), e.vertex(1), e.vertex(2) });
										}

										VolumeMesh<double> *vm = new VolumeMesh<double>(std::move(volume_mesh_elements),
										                                                std::move(volume_mesh_vertices));
										VolumeMeshFieldLinear<double, double> *pf = new VolumeMeshFieldLinear<double, double>(
										    MakeConvexPressureField<double>(vm, hydroElasticModulus));
										Bvh<Obb, VolumeMesh<double>> *bvh = new Bvh<Obb, VolumeMesh<double>>(*vm);
										cp = new ContactProperties(id, s, SOFT, NULL, vm, pf, bvh, hydroElasticModulus,
										                           dissipation, staticFriction, dynamicFriction);
									} else {
										Bvh<Obb, TriangleSurfaceMesh<double>> *bvh =
										    new Bvh<Obb, TriangleSurfaceMesh<double>>(*sm);
										cp = new ContactProperties(id, s, RIGID, NULL, sm, bvh, staticFriction, dynamicFriction);
									}
									contactProperties[id] = cp;
									break;
								}
							}
							ROS_WARN_STREAM_NAMED("mujoco_contact_surfaces",
							                      "Could not load mesh! It was not properly defined in mujoco.");
							break;
						}
					}
				}
			}
		}
	}
}

void MujocoContactSurfacesPlugin::renderCallback(mjModelPtr model, mjDataPtr data, mjvScene *scene)
{
	if (visualizeContactSurfaces) {
		int n = std::min(n_vGeom, scene->maxgeom);
		for (int i = 0; i < n; ++i) {
			scene->geoms[scene->ngeom++] = vGeoms[i];
		}
		n_vGeom = 0;
	}
}

} // namespace mujoco_contact_surfaces

PLUGINLIB_EXPORT_CLASS(mujoco_contact_surfaces::MujocoContactSurfacesPlugin, MujocoSim::MujocoPlugin)
