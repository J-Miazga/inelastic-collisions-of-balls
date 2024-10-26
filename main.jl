# using GeometryBasics
# using GLMakie
# using DifferentialEquations
# using OrdinaryDiffEq

# Function to generate the coordinates for a circle

# mutable struct Ball
#     position::Vector{Float64}
#     force::Vector{Float64}
#     mass::Float64
#     radius::Float64
# end

# using GLMakie
# using Observables
# using GeometryBasics

# mutable struct Ball
#     mass::Float32
#     position::Vector{Float32}
#     force::Vector{Float32}
# end
# function generate_circle_points(center, radius::Float64, n_points = 100)
#     angles = LinRange(0, 2pi, n_points)  # Create angle steps
#     xs = center[1] .+ radius .* cos.(angles)  # X coordinates
#     ys = center[2] .+ radius .* sin.(angles)  # Y coordinates
#     return xs, ys
# end
# balls=Vector{Ball}()
# # Create a figure and axis
# fig = Figure(size=(1920, 1080))
# ax = Axis(fig[1, 1],aspect=DataAspect())
# Makie.deactivate_interaction!(ax,:rectanglezoom)

# # Create a line plot

# lines!(ax,Circle(Point2f0(1.0,1.0), 1.0), color = :black, linewidth = 10)
# buttons_grid= fig[2,1]= GridLayout(tellwidth=false)
# return mass_box = Textbox(buttons_grid[1,1],placeholder="Mass", height=30,width=150)
# force_box = Textbox(buttons_grid[1,2],placeholder="Force",height=30,width=150)
# button1 = Button(buttons_grid[1,3], label="Submit",height=30,width=150)
# Label(buttons_grid[1, 4], "Run Simulation")
# toggle=Toggle(buttons_grid[1,5],height=30,width=50)

# spoint =select_point(ax.scene)
# control=1
# position_of_center=Vector{Float32}()
# alfa=0
# on(spoint) do position
#     if control%2==1
#         x,y=generate_circle_points(position, 0.1, 100)
#         poly!(ax,x,y, color = :black)
#         println("Position x= $(position[1]) and y= $(position[2])")
#         global position_of_center=position
#         global control+=1
#     else
#         global alfa=atand(abs(position[2]-position_of_center[2])/abs(position[1]-position_of_center[1]))
#     end
# end 
# on(button1.clicks) do g
#     if !isnothing(mass_box.stored_string[]) && !isnothing(force_box.stored_string[])
#         vector_force=Vector{Float32}()
#         mass = tryparse(Float32, mass_box.stored_string[])
#         force = tryparse(Float32, force_box.stored_string[])
#         vector_force=[force*cos(alfa),force*sin(alfa)]
#         new_ball=Ball(mass,position_of_center,vector_force)
#         push!(balls,new_ball)
#         global control+=1
#         mass_box.stored_string=nothing
#         force_box.stored_string=nothing
#     end
# end

# # Display the plot
# display(fig)
using GLMakie
using Observables
using GeometryBasics

mutable struct Ball
    mass::Float32
    position::Vector{Float32}
    force::Vector{Float32}
end
mutable struct Balls_State
    control::Int
    position_of_center::Vector{Float32}
    alfa::Float32
end
function init_state()
    init=Balls_State(1,[0.0,0.0],0.0)
end
function generate_circle(center, radius::Float64, n_points = 100)
    angles = LinRange(0, 2pi, n_points) 
    xs = center[1] .+ radius .* cos.(angles)  
    ys = center[2] .+ radius .* sin.(angles)  
    return xs, ys
end
function click_plot_check(state::Balls_State, ax, position)
    if state.control % 2 == 1
        x, y = generate_circle(position, 0.1, 100)
        poly!(ax, x, y, color=:black)
        state.position_of_center = position  
        state.control += 1                 
    else
        state.alfa = atand(abs(position[2] - state.position_of_center[2]) / abs(position[1] - state.position_of_center[1]))
    end
end
function button_click_check(state::Balls_State, balls, mass_box, force_box)
    if !isnothing(mass_box.stored_string[]) && !isnothing(force_box.stored_string[])
        mass = tryparse(Float32, mass_box.stored_string[])
        force = tryparse(Float32, force_box.stored_string[])
        vector_force = [force * cos(state.alfa), force * sin(state.alfa)]       
        new_ball = Ball(mass, state.position_of_center, vector_force)
        push!(balls, new_ball)

        state.control += 1
        mass_box.stored_string= nothing
        force_box.stored_string= nothing
    end
end

create_state = init_state()
fig = Figure(size=(1920, 1080))
ax = Axis(fig[1, 1],aspect=DataAspect())
Makie.deactivate_interaction!(ax,:rectanglezoom)

lines!(ax,Circle(Point2f0(1.0,1.0), 1.0), color = :black, linewidth = 10)
buttons_grid= fig[2,1]= GridLayout(tellwidth=false)
mass_box = Textbox(buttons_grid[1,1],placeholder="Mass",reset_on_defocus=true, height=30,width=150)
force_box = Textbox(buttons_grid[1,2],placeholder="Force",reset_on_defocus=true,height=30,width=150)
button1 = Button(buttons_grid[1,3], label="Submit",height=30,width=150)
Label(buttons_grid[1, 4], "Run Simulation")
toggle=Toggle(buttons_grid[1,5],height=30,width=50)
spoint = select_point(ax.scene)
balls=Vector{Ball}()

on(spoint) do position
    click_plot_check(create_state, ax, position)
end

on(button1.clicks) do g
    button_click_check(create_state, balls, mass_box, force_box)
end

display(fig)

























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
