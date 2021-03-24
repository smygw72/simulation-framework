#pragma warning( disable : 4996)

#include "simulation.h"
#include "timer_wrapper.h"

Simulation::Simulation()
{
}

Simulation::~Simulation()
{
	clearConstraints();
}

void Simulation::reset()
{    
	std::cout << "Reset Simulation ..." << std::endl;

	// set pointer of m_V, m_F, m_T
	m_V = &m_mesh->m_V;
	m_F = &m_mesh->m_F;
	m_T = &m_mesh->m_T;
	m_Vel.resize(m_mesh->m_vert_num, 3); m_Vel.setZero();

	m_Inertia.resize(m_mesh->m_vert_num, 3);       m_Inertia.setZero();
	m_ExternalForce.resize(m_mesh->m_vert_num, 3); m_ExternalForce.setZero();

	setReprecomputeFlag();
	setReprefactorFlag();

	setupConstraints();
	preComputation();

	m_selected_attachment_constraint = NULL;
}

void Simulation::update()
{
	// update inertia term
	computeInertia();

	// update external force
	computeExternalForce();

	// update cloth
	switch (m_integration_method)
	{
	case INTEGRATION_LOCAL_GLOBAL:
		integrateOptimizationMethod();
		break;
	}

	// Add collision detection here using pos_next;
	// Inprement Here !! //
	// Hint: use "collisionDetection"
	// This is optional (you don't have to implement first)

	// update velocity and damp
	dampVelocity();
}

void Simulation::convertLameConstant()
{
	m_myu = m_young / (2.f * (1.f + m_poisson));
	m_lambda = m_young * m_poisson / ((1.f + m_poisson)*(1.f - 2.f*m_poisson));
}

#pragma region Constraint
void Simulation::drawConstraints(const VBO& vbos)
{
	for (std::vector<Constraint*>::iterator it = m_constraints.begin(); it != m_constraints.end(); ++it)
	{
		(*it)->draw(vbos);
	}
}

ScalarType Simulation::tryToSelectAttachmentConstraint(const EigenVector3& p0, const EigenVector3& dir)
{
	ScalarType ray_point_dist;
	ScalarType min_dist = 100.0;
	AttachmentConstraint* best_candidate = NULL;

	bool current_state_on = false;
	for (std::vector<Constraint*>::iterator c = m_constraints.begin(); c != m_constraints.end(); ++c)
	{
		AttachmentConstraint* ac;
		if (ac = dynamic_cast<AttachmentConstraint*>(*c)) // is attachment constraint
		{
			ray_point_dist = ((ac->getFixedPoint()-p0).cross(dir)).norm();
			if (ray_point_dist < min_dist)
			{
				min_dist = ray_point_dist;
				best_candidate = ac;
			}
		}
	}
	// exit if no one fits
	if (min_dist > DEFAULT_SELECTION_RADIUS)
	{
		unselectAttachmentConstraint();

		return -1;
	}
	else
	{
		selectAtttachmentConstraint(best_candidate);
		EigenVector3 fixed_point_temp = m_V->row(m_selected_attachment_constraint->getConstrainedVertexIndex()).transpose();

		return (fixed_point_temp-p0).dot(dir); // this is m_cached_projection_plane_distance
	}
}

bool Simulation::tryToToggleAttachmentConstraint(const EigenVector3& p0, const EigenVector3& dir)
{
	EigenVector3 p1;

	ScalarType ray_point_dist;
	ScalarType min_dist = 100.0;
	unsigned int best_candidate = 0;
	// first pass: choose nearest point
	for (unsigned int i = 0; i != m_mesh->m_vert_num; i++)
	{
		p1 = m_V->row(i).transpose();
		
		ray_point_dist = ((p1-p0).cross(dir)).norm();
		if (ray_point_dist < min_dist)
		{
			min_dist = ray_point_dist;
			best_candidate = i;
		}
	}
	for (std::vector<Constraint*>::iterator c = m_constraints.begin(); c != m_constraints.end(); ++c)
	{
		AttachmentConstraint* ac;
		if (ac = dynamic_cast<AttachmentConstraint*>(*c)) // is attachment constraint
		{
			ray_point_dist = ((ac->getFixedPoint()-p0).cross(dir)).norm();
			if (ray_point_dist < min_dist)
			{
				min_dist = ray_point_dist;
				best_candidate = ac->getConstrainedVertexIndex();
			}
		}
	}
	// exit if no one fits
	if (min_dist > DEFAULT_SELECTION_RADIUS)
	{
		return false;
	}
	// second pass: toggle that point's fixed position constraint
	bool current_state_on = false;
	for (std::vector<Constraint*>::iterator c = m_constraints.begin(); c != m_constraints.end(); ++c)
	{
		AttachmentConstraint* ac;
		if (ac = dynamic_cast<AttachmentConstraint*>(*c)) // is attachment constraint
		{
			if (ac->getConstrainedVertexIndex() == best_candidate)
			{
				current_state_on = true;
				m_constraints.erase(c);
				break;
			}
		}
	}
	if (!current_state_on)
	{
		addAttachmentConstraint(best_candidate);
	}

	return true;
}

void Simulation::selectAtttachmentConstraint(AttachmentConstraint* ac)
{
	m_selected_attachment_constraint = ac;
	m_selected_attachment_constraint->select();
}

void Simulation::unselectAttachmentConstraint()
{
	if (m_selected_attachment_constraint)
	{
		m_selected_attachment_constraint->unSelect();
	}
	m_selected_attachment_constraint = NULL;
}

void Simulation::addAttachmentConstraint(uint vertex_index)
{
	EigenVector3 p = m_mesh->m_V.row(vertex_index).transpose();
	AttachmentConstraint* ac = new AttachmentConstraint(&m_stiffness_attachment, vertex_index, p);
	m_constraints.push_back(ac);
}

void Simulation::moveSelectedAttachmentConstraintTo(const EigenVector3& target)
{
	if (m_selected_attachment_constraint)
		m_selected_attachment_constraint->setFixedPoint(target);
}

void Simulation::clearConstraints()
{
	for (unsigned int i = 0; i < m_constraints.size(); ++i)
	{
		delete m_constraints[i];
	}
	m_constraints.clear();
}

void Simulation::setupConstraints()
{
	clearConstraints();

	// set init constraint
	switch(m_mesh->m_mesh_type)
	{
	case MESH_TYPE_TET:
		{
			// set no constraint
		}
		break;
	}
}
#pragma endregion

#include<igl/massmatrix.h>
void Simulation::preComputation()
{
	std::cout << "preComputing..." << std::endl;

	// set m_B, m_W
	m_B.clear();
	m_W.clear();

	// Inprement Here !! //
	// Hint: compute and push_back "m_B" and "m_W"
	// 

	// set MassMatrix
	/*igl::massmatrix(*m_V, *m_T, igl::MASSMATRIX_TYPE_DEFAULT, m_MassMat);
	ScalarType ModelVolume ;
	for (uint i = 0; i<m_T->rows(); ++i) { ModelVolume += m_W[i]; }
	m_MassMat *= (m_mesh->m_total_mass) / ModelVolume;*/

	m_precomputation_flag = true;
}

void Simulation::dampVelocity()
{
	if (std::abs(m_damping_coefficient) < EPSILON)
		return;

	// post-processing damping
	EigenVector3 pos_mc(0.0, 0.0, 0.0), vel_mc(0.0, 0.0, 0.0);
	unsigned int i, size;
	ScalarType denominator(0.0), mass(0.0);
	size = m_mesh->m_vert_num;

	{
		for (int i = 0; i < size; ++i)
		{
			mass = m_MassMat.coeff(i, i);

			pos_mc += mass*(((*m_V).row(i)).transpose());
			vel_mc += mass*(((m_Vel).row(i)).transpose());
			denominator += mass;
		}

		assert(denominator != 0.0);
		pos_mc /= denominator;
		vel_mc /= denominator;

		EigenVector3 angular_momentum(0.0, 0.0, 0.0), r(0.0, 0.0, 0.0);
		EigenMatrix3 inertia, r_mat;
		inertia.setZero(); r_mat.setZero();

		for (int i = 0; i < size; ++i)
		{
			mass = m_MassMat.coeff(i, i);

			r = ((*m_V).row(i)).transpose() - pos_mc;


			//r_mat = EigenMatrix3(0.0,  r.z, -r.y,
			//                    -r.z, 0.0,  r.x,
			//                    r.y, -r.x, 0.0);

			r_mat.coeffRef(0, 1) = r[2];
			r_mat.coeffRef(0, 2) = -r[1];
			r_mat.coeffRef(1, 0) = -r[2];
			r_mat.coeffRef(1, 2) = r[0];
			r_mat.coeffRef(2, 0) = r[1];
			r_mat.coeffRef(2, 1) = -r[0];

			inertia += r_mat * r_mat.transpose() * mass;
		}
		EigenVector3 angular_vel = inertia.inverse() * angular_momentum;

		EigenVector3 delta_v(0.0, 0.0, 0.0);

		for (int i = 0; i < size; ++i)
		{
			r = ((*m_V).row(i)).transpose() - pos_mc;
			delta_v = vel_mc + angular_vel.cross(r) - ((m_Vel.row(i)).transpose());
			m_Vel.row(i) += m_damping_coefficient * (delta_v.transpose());
		}
	}
}

void Simulation::computeInertia()
{
	// Inprement Here !! //
	// Hint: use "m_Inertia = ?"
}

void Simulation::computeExternalForce()
{
	
	m_ExternalForce.setZero();

	// gravity
	for (unsigned int i = 0; i < m_mesh->m_vert_num; ++i)
	{
		// Inprement Here !! //
		// Hint: use "m_ExternalForce(i, 1) = ?"
	}
}

void Simulation::collisionDetection(EigenMatrixXs& x, EigenMatrixXs& v)
{
	// Inprement Here !! //
	// optional part!
}

#pragma region integration scheme
void Simulation::integrateOptimizationMethod()
{
	// check if precomputation is done (for Debag)
	if (m_precomputation_flag == false) { fprintf(stdout, "precompution dosen't work\n"); }

	// take a initial guess
	EigenMatrixXs pos_next = m_Inertia;

	// while loop until converge or exceeds maximum iterations
	bool converge = false;

	for (unsigned int iteration_num = 0; !converge && iteration_num < m_iterations_per_frame; ++iteration_num)
	{
		switch (m_integration_method)
		{
		case INTEGRATION_LOCAL_GLOBAL:
			converge = integrateLocalGlobalOneIteration(pos_next);
			break;
		}
	}

	// update q_{n+1}
	updatePosAndVel(pos_next);
}

bool Simulation::integrateLocalGlobalOneIteration(EigenMatrixXs& X)
{
	prefactorize();

	// local step
	EigenMatrixXs RotMat(m_T->rows() * 3, 3);       // (#TetNum*dim) * dim
	EigenMatrixXs Jv = m_JacobianMat * X;      // (#TetNum*dim) * dim
	computeRotMat(RotMat, Jv);

	// global step
	EigenMatrixXs b(m_mesh->m_vert_num, 3);
	// Inprement Here !! //
	// Hint: use "b = ?"
	// Hint: you should not add inertia yet

	// add attachment constraint
	for (std::vector<Constraint*>::iterator it_c = m_constraints.begin(); it_c != m_constraints.end(); ++it_c)
	{
		(*it_c)->computeJVector(X, b);
	}

	// add effect of inertia
	// Inprement Here !! //
	// Hint: use "b = ?"
	// Hint: add intertia here

	X = m_prefactored_LLTsolver.solve(b);

	return false;
}
#pragma endregion

#pragma region local_global
void Simulation::computeRotMat(EigenMatrixXs& RotMat, const EigenMatrixXs& Jv)
{
	// Inprement Here !! //
	// Hint: use "SVD"
}

void Simulation::computeElementLaplacianMat(const EigenMatrix3 &B, const ScalarType W, const unsigned int tet_list[], std::vector<SparseMatrixTriplet>& l_triplets)
{
	// Inprement Here !! //
	// Small Hint: use "m_myu"
}

void Simulation::computeElementJacobianMat(const EigenMatrix3 &B, const ScalarType W, const unsigned int tet_list[], const unsigned int ele_num, std::vector<SparseMatrixTriplet>& j_triplets)
{
	// Inprement Here !! //
	// Small Hint: use "m_myu"
}
#pragma endregion

#pragma region matrices and prefactorization
void Simulation::setLaplacianMat()
{
	m_LaplacianMat.resize(m_mesh->m_system_dim, m_mesh->m_system_dim);
	std::vector<SparseMatrixTriplet> l_triplets;
	l_triplets.clear();

	for (uint i = 0; i<m_T->rows(); ++i)
	{
		uint tet_list[4] = { (*m_T)(i, 0), (*m_T)(i, 1),(*m_T)(i, 2), (*m_T)(i, 3) };
		// Inprement Here !! //
		// Hint: use "computeElementLaplacianMat"
	}

	// add attachment constraint
	for (std::vector<Constraint*>::iterator it_c = m_constraints.begin(); it_c != m_constraints.end(); ++it_c)
	{
		(*it_c)->computeLaplacianMat(l_triplets);
	}

	m_LaplacianMat.setFromTriplets(l_triplets.begin(), l_triplets.end());
}

void Simulation::setJacobianMat()
{
	m_JacobianMat.resize(m_T->rows() * 3, m_mesh->m_vert_num);
	std::vector<SparseMatrixTriplet> j_triplets;
	j_triplets.clear();

	for (uint i = 0; i<m_T->rows(); ++i)
	{
		uint tet_list[4] = { (*m_T)(i,0), (*m_T)(i,1), (*m_T)(i,2), (*m_T)(i,3) };
		// Inprement Here !! //
		// Hint: use "computeElementJacobianMat"
	}

	m_JacobianMat.setFromTriplets(j_triplets.begin(), j_triplets.end());
}

void Simulation::prefactorize()
{
	SparseMatrix A;
	ScalarType h2 = m_h*m_h;

	std::cout << "preFactorizing ..." << std::endl;

	// Inprement Here !! //
	// Hint: "A = ?"

	//factorizeDirectSolverLLT(A, m_prefactored_LLTsolver);

	m_prefactorization_flag = true;
}
#pragma endregion

#pragma region utilities
void Simulation::updatePosAndVel(const EigenMatrixXs& NewPos)
{
	m_Vel = (NewPos - (*m_V)) / m_h;
	*m_V  = NewPos;
}
#pragma endregion
