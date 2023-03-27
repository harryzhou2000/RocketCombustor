PROGRAM = RocketCombustor
$(PROGRAM): module.f90 \
			parallel.f90 \
			combustion_modeling/species.f90 \
			combustion_modeling/fuel.f90 \
			combustion_modeling/oxidizer.f90 \
			combustion_modeling/thermo_property.f90 \
			combustion_modeling/transport_property.f90 \
			combustion_modeling/diffusion_velocity.f90 \
			combustion_modeling/chemical_kinetics.f90 \
			rnn.f90 \
			period.f90 \
			allocate_arrary.f90 \
			flux_lax.f90 \
			flux_viscous.f90 \
			fvm2d.f90 \
			gaussInteg4.f90 \
			ghostcell.f90 \
			gradient.f90 \
			mesh.f90 \
			initialization.f90 \
			initParam.f90 \
			output.f90 \
			reconstruction_initialization.f90 \
			limiter.f90 \
			shockdetect.f90 \
			startup.f90 \
			Time_RK_LUSGS.f90 \
			Time_LUSGS_steady.f90 \
			Time_RungeKutta.f90 \
			TimeStep.f90 \
			sourceterm.f90 \
			residual.f90 \
			monitor.f90 \
			common_boundary.f90 \
			update_flow_variables.f90 \
			la_math.f90 \
			Con_to_Pri_Jacobian.f90 \
			Time_Trapezoidal.f90 \
			Time_BackwardEuler.f90 \
			Time_ImplicitMidpoint.f90 \
			Time_BDF2.f90
	mpifort -o $@ $^ -O3 -lblas -llapack -fdefault-double-8 -fdefault-real-8 
# mpiifort -lblas -llapcak -O3 -autodouble -fpconstant -o $@ $^
# mpifort -o $@ $^ -O3 -lblas -llapack -fdefault-double-8 -fdefault-real-8 
	

.PHONY:clean
clean:
	rm -rf $(PROGRAM) *.mod *.o *.vtk *.dat *.plt *.visit
