PROGRAM = meshpre
$(PROGRAM): module.f90 \
			period.f90 \
			allocate_arrary.f90 \
			pre-process.f90 \
			geometry_mix.f90 \
			MeshPre.f90 \
			metisgrid.f90 \
			reconst_init.f90

	gfortran -fcheck=all -fdefault-double-8 -fdefault-real-8 -o $@ $^
	
.PHONY:clean
clean:
	rm -rf $(PROGRAM) *.mod *.o *.vtk *.dat *.plt *.visit *.PLT
