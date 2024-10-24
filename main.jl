# using GeometryBasics
# using GLMakie
# using DifferentialEquations
# using OrdinaryDiffEq

# Function to generate the coordinates for a circle
# function generate_circle_points(center, radius, n_points = 100)
#     angles = LinRange(0, 2pi, n_points)  # Create angle steps
#     xs = center[1] .+ radius .* cos.(angles)  # X coordinates
#     ys = center[2] .+ radius .* sin.(angles)  # Y coordinates
#     return xs, ys
# end
# mutable struct Ball
#     position::Vector{Float64}
#     force::Vector{Float64}
#     mass::Float64
#     radius::Float64
# end

using GLMakie
using Observables
using GeometryBasics

# Create a figure and axis
fig = Figure(size=(1920, 1080))
ax = Axis(fig[1, 1],aspect=DataAspect())
Makie.deactivate_interaction!(ax,:rectanglezoom)
# Create a line plot

lines!(ax,Circle(Point2f0(1.0,1.0), 1.0), color = :black, linewidth = 10)
number_of_balls=Observable(0)
# Adding a textbox for number of balls
text_box = Textbox(fig[2,1],placeholder="Number of balls", tellwidth = false)

# Create a button to confirm the input
button1 = Button(fig[2,2], label="Submit")
button2 = Button(fig[2,3], label="Run")

on(button1.clicks) do n
    println("Button clicked")
    number_of_balls[] = tryparse(Int, text_box.stored_string[])
    println("Number of balls set to: $(number_of_balls[])")  # Show valid integer   
end
ball_mat=Matrix{Float32}(undef,0,2)
spoint =select_point(ax.scene)
on(spoint) do z
    if number_of_balls[]!=size(ball_mat,1)
    global ball_mat
    position=[0.0,0.0]
    position=z
    ball_mat=vcat(ball_mat,position')
    lines!(ax,Circle(Point2f0(ball_mat[size(ball_mat,1),1],ball_mat[size(ball_mat,1),2]), 0.1), color = :black,linewidth = 2)
    println("Position x= $(position[1]) and y= $(position[2])")
    end
end

# Display the plot
display(fig)
 
# fig = Figure(size = (1920, 1080))
# ax = Axis(fig[1, 1], aspect= DataAspect())
# Makie.deactivate_interaction!(ax,:rectanglezoom)
# spoint =select_point(ax.scene)


# on(spoint) do z
#     position_x,positon_y=z
#     println(position_x)
# end

# display(fig)

























# matrix_of_balls=[]
# push!(matrix_of_balls,[10.0,1.0,20.0,20.0,4.0,0.1])


# v_initial=[ball1.force[1]/ball1.mass,ball1.force[2]/ball1.mass]
# v_x_initial=v_initial[1].*cos(rad2deg(atan(matrix_of_balls[1,1]/matrix_of_balls[1,2])))
# v_y_initial=v_initial[2].*sin(rad2deg(atan(ball1.position[1]/ball1.position[2])))
# u0_matrix=[]
# push!(u0_matrix,[])
# u0 = [v_x_initial,v_y_initial, ball1.position[1],ball1.position[2]]
# p=[c,ball1.mass,g]
# function ball_motion!(du,u,p,t)
#     v_x, v_y, x,y=u
#     c,m,g=p
#     du[1] = - c/m *  v_x  # dvx/dt (air resistance in x)
#     du[2] = g - c/m * v_y  # dvy/dt (gravity + air resistance)
#     du[3] = v_x  # dx/dt = vx
#     du[4] = v_y  # dy/dt = vy
   
# end
# integ=integrator(ball_motion!,u0,Tsit5(),dt=0.005)
# fig = Figure(size = (1920, 1080))
# ax = Axis(fig[1, 1], aspect= DataAspect())

# tspan = (0.0, 10.0)
# prob = ODEProblem(ball_motion!, u0, tspan,p)
# sol=solve(prob,Tsit5();saveat=0.01)

# # Define the circle parameters
# #center = Point2f0(0.0, 0.0)  # Circle center
# #radius = 1.0  # Circle radius
# # Generate circle points
# #xs, ys = generate_circle_points(center, radius)

# # Plot the circle using lines
# #lines!(ax, xs, ys, color = :black, linewidth = ox)
# #lines!(ax,Circle(Point2f0(0.0,0.0), 1.0), color = :black, linewidth = 10)
# lines!(ax,Point2f0(ball1.position), color = :black)

# # Show the figure
# display(fig)
# #empty!(ax)
# # while events(fig).window_open[
# #     display(fig)
# # end
