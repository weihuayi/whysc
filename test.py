import os

<<<<<<< HEAD
a = [1, 4, 16]
for i in range(3):
    os.system("build/example/pmesh_generation quad_model1 %i"%(a[i]))
    os.system("mpirun -n %i build/example/parallel_quadmesh_optimization test_quad"%(a[i]))

    os.system("build/example/pmesh_trimesh_generation tri_model1 %i"%(a[i]))
    os.system("mpirun -n %i build/example/parallel_trimesh_optimization test_tri"%(a[i]))

    os.system("build/example/pmesh_tetmesh_generation tet_model1 %i"%(a[i]))
    os.system("mpirun -n %i build/example/parallel_tetmesh_optimization test_tet"%(a[i]))

    os.system("build/example/pmesh_hexmesh_generation hex_model1 %i"%(a[i]))
    os.system("mpirun -n %i build/example/parallel_hexmesh_optimization test_hex"%(a[i]))


    os.system("build/example/pmesh_generation quad_model2 %i"%(a[i]))
    os.system("mpirun -n %i build/example/parallel_quadmesh_optimization0 test_quad"%(a[i]))

    os.system("build/example/pmesh_trimesh_generation tri_model2 %i"%(a[i]))
    os.system("mpirun -n %i build/example/parallel_trimesh_optimization0 test_tri"%(a[i]))

    os.system("build/example/pmesh_tetmesh_generation tet_model2 %i"%(a[i]))
    os.system("mpirun -n %i build/example/parallel_tetmesh_optimization0 test_tet"%(a[i]))

    os.system("build/example/pmesh_hexmesh_generation hex_model2 %i"%(a[i]))
=======
a = np.array([1, 4, 16])
for i in range(3):
    os.system("./build/pmesh_generation quad_model1 %i"%(a[i]))
    os.system("mpirun -n %i build/example/parallel_quadmesh_optimization test_quad"%(a[i]))

    os.system("./build/pmesh_trimesh_generation tri_model1 %i"%(a[i]))
    os.system("mpirun -n %i build/example/parallel_trimesh_optimization test_tri"%(a[i]))

    os.system("./build/pmesh_tetmesh_generation tet_model1 %i"%(a[i]))
    os.system("mpirun -n %i build/example/parallel_tetmesh_optimization test_tet"%(a[i]))

    os.system("./build/pmesh_hexmesh_generation hex_model1 %i"%(a[i]))
    os.system("mpirun -n %i build/example/parallel_hexmesh_optimization test_hex"%(a[i]))


    os.system("./build/pmesh_generation quad_model2 %i"%(a[i]))
    os.system("mpirun -n %i build/example/parallel_quadmesh_optimization0 test_quad"%(a[i]))

    os.system("./build/pmesh_trimesh_generation tri_model2 %i"%(a[i]))
    os.system("mpirun -n %i build/example/parallel_trimesh_optimization0 test_tri"%(a[i]))

    os.system("./build/pmesh_tetmesh_generation tet_model2 %i"%(a[i]))
    os.system("mpirun -n %i build/example/parallel_tetmesh_optimization0 test_tet"%(a[i]))

    os.system("./build/pmesh_hexmesh_generation hex_model2 %i"%(a[i]))
>>>>>>> a4b3490af0ba760e007941c391ed5fc820ad9433
    os.system("mpirun -n %i build/example/parallel_hexmesh_optimization0 test_hex"%(a[i]))
