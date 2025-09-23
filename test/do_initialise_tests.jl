
using TOML, Skyrmions3D

function test_if_skyrmions_equal(a_skyrmion, b_skyrmion)

    @test a_skyrmion.pion_field == b_skyrmion.pion_field
    @test a_skyrmion.grid.lp == b_skyrmion.grid.lp
    @test a_skyrmion.grid.ls == b_skyrmion.grid.ls
    @test a_skyrmion.mpi == b_skyrmion.mpi
    @test a_skyrmion.Fpi == b_skyrmion.Fpi
    @test a_skyrmion.ee == b_skyrmion.ee
    @test a_skyrmion.physical == b_skyrmion.physical
    @test a_skyrmion.grid.dirichlet == b_skyrmion.grid.dirichlet
    @test a_skyrmion.grid.index_grid_x == b_skyrmion.grid.index_grid_x
    @test a_skyrmion.grid.index_grid_y == b_skyrmion.grid.index_grid_y
    @test a_skyrmion.grid.index_grid_z == b_skyrmion.grid.index_grid_z
    @test a_skyrmion.grid.sum_grid == b_skyrmion.grid.sum_grid
    @test a_skyrmion.vac == b_skyrmion.vac

end

a_skyrmion = Skyrmion(5, 0.2)
the_same_skyrmion = Skyrmion([5, 5, 5], [0.2, 0.2, 0.2])

test_if_skyrmions_equal(a_skyrmion, the_same_skyrmion)

# setting stuff

set_mpi!(a_skyrmion, 1.0)
@test a_skyrmion.mpi == 1.0

set_Fpi!(a_skyrmion, 100)
@test a_skyrmion.Fpi == 100

set_ee!(a_skyrmion, 3.0)
@test a_skyrmion.ee == 3.0

set_lattice!(a_skyrmion, [6, 4, 7], [0.1, 0.15, 0.25])

@test a_skyrmion.grid.lp[1] == 6
@test a_skyrmion.grid.lp[2] == 4
@test a_skyrmion.grid.lp[3] == 7

@test a_skyrmion.grid.ls[1] == 0.1
@test a_skyrmion.grid.ls[2] == 0.15
@test a_skyrmion.grid.ls[3] == 0.25

@test a_skyrmion.grid.dirichlet == true

set_neumann!(a_skyrmion)
@test a_skyrmion.grid.dirichlet == false
@test a_skyrmion.grid.boundary_conditions == "neumann"

@test a_skyrmion.grid.sum_grid == [1:6, 1:4, 1:7]

@test a_skyrmion.grid.index_grid_x == [2, 1, 1, 2, 3, 4, 5, 6, 6, 5]
@test a_skyrmion.grid.index_grid_y == [2, 1, 1, 2, 3, 4, 4, 3]
@test a_skyrmion.grid.index_grid_z == [2, 1, 1, 2, 3, 4, 5, 6, 7, 7, 6]

set_periodic!(a_skyrmion)
@test a_skyrmion.grid.dirichlet == false
@test a_skyrmion.grid.boundary_conditions == "periodic"

@test a_skyrmion.grid.sum_grid == [1:6, 1:4, 1:7]

@test a_skyrmion.grid.index_grid_x == [5, 6, 1, 2, 3, 4, 5, 6, 1, 2]
@test a_skyrmion.grid.index_grid_y == [3, 4, 1, 2, 3, 4, 1, 2]
@test a_skyrmion.grid.index_grid_z == [6, 7, 1, 2, 3, 4, 5, 6, 7, 1, 2]

set_dirichlet!(a_skyrmion)
@test a_skyrmion.grid.dirichlet == true
@test a_skyrmion.grid.boundary_conditions == "dirichlet"

@test a_skyrmion.grid.sum_grid == [3:4, 3:2, 3:5]

set_physical!(a_skyrmion, true)

@test a_skyrmion.physical == true

set_physical!(a_skyrmion, false)

@test a_skyrmion.physical == false

@test Skyrmions3D.sum_grid([3, 3, 3], "dirichlet") == Skyrmions3D.sum_grid(3, "dirichlet")

@test Skyrmions3D.sum_grid([3, 3, 3], "neumann") == Skyrmions3D.sum_grid(3, "neumann")

@test Skyrmions3D.sum_grid([3, 3, 3], "periodic") == Skyrmions3D.sum_grid(3, "periodic")

@test_logs (:warn, "Unrecognised boundary conditions: unexpected behaviour may occur") Skyrmions3D.sum_grid(
    [3, 3, 3],
    "nonsense_boundary_condition",
)

@test Skyrmions3D.index_grid(6, "periodic")[end] == 2
@test Skyrmions3D.index_grid(6, "periodic")[end-1] == 1
@test Skyrmions3D.index_grid(6, "periodic")[1] == 5
@test Skyrmions3D.index_grid(6, "periodic")[2] == 6

a_skyrmion.pion_field[2, 3, 1, 1] = 2.0

@test_throws AssertionError check_if_normalised(a_skyrmion)

normer!(a_skyrmion)
check_if_normalised(a_skyrmion)

a_skyrmion.pion_field[2, 3, 1, 1] = 2.0
check_if_normalised(normer(a_skyrmion))


# wrapper to make a temp directory, which will deleted when tests are completed.
mktempdir() do tmpdir

    filepath = joinpath(tmpdir, "basic_save")
    save_skyrmion(a_skyrmion, filepath)
    loaded_skyrmion = load_skyrmion(filepath)
    test_if_skyrmions_equal(a_skyrmion, loaded_skyrmion)

    filepath = joinpath(tmpdir, "with_user_metadata")
    save_skyrmion(a_skyrmion, filepath, additional_metadata = Dict("some" => "metadata"))
    loaded_skyrmion = load_skyrmion(filepath)
    test_if_skyrmions_equal(a_skyrmion, loaded_skyrmion)

    metadata = nothing
    open(joinpath(filepath, "metadata.toml"), "r") do io
        metadata = TOML.parse(io)
    end
    @test metadata["additional_metadata"]["some"] == "metadata"


end
