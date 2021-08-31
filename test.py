import os

a = [4, 8]
for i in range(2):
    print("第 %i 步, 第 1 次"%i)
    os.system("build/example/pmesh_generation quad_model1 %i"%(a[i]))
    os.system("mpirun -n %i build/example/parallel_quadmesh_optimization test_quad"%(a[i]))

    print("第 %i 步, 第 2 次"%i)
    os.system("build/example/pmesh_trimesh_generation tri_model1 %i"%(a[i]))
    os.system("mpirun -n %i build/example/parallel_trimesh_optimization test_tri"%(a[i]))

    print("第 %i 步, 第 3 次"%i)
    os.system("build/example/pmesh_tetmesh_generation tet_model1 %i"%(a[i]))
    os.system("mpirun -n %i build/example/parallel_tetmesh_optimization test_tet"%(a[i]))

    print("第 %i 步, 第 4 次"%i)
    os.system("build/example/pmesh_hexmesh_generation hex_model1 %i"%(a[i]))
    os.system("mpirun -n %i build/example/parallel_hexmesh_optimization test_hex"%(a[i]))


    print("第 %i 步, 第 5 次"%i)
    os.system("build/example/pmesh_generation quad_model2 %i"%(a[i]))
    os.system("mpirun -n %i build/example/parallel_quadmesh_optimization0 test_quad"%(a[i]))

    print("第 %i 步, 第 6 次"%i)
    os.system("build/example/pmesh_trimesh_generation tri_model2 %i"%(a[i]))
    os.system("mpirun -n %i build/example/parallel_trimesh_optimization0 test_tri"%(a[i]))

    print("第 %i 步, 第 7 次"%i)
    os.system("build/example/pmesh_tetmesh_generation tet_model2 %i"%(a[i]))
    os.system("mpirun -n %i build/example/parallel_tetmesh_optimization0 test_tet"%(a[i]))

    print("第 %i 步, 第 8 次"%i)
    os.system("build/example/pmesh_hexmesh_generation hex_model2 %i"%(a[i]))
    os.system("mpirun -n %i build/example/parallel_hexmesh_optimization0 test_hex"%(a[i]))
