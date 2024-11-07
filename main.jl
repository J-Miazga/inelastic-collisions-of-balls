using GLMakie
using Observables
using GeometryBasics
using LinearAlgebra
using Dates

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
function generate_circle(center::Vector{Float32}, radius::Float64, n_points::Int)
    angles = LinRange(0, 2pi, n_points) 
    xs = center[1] .+ radius .* cos.(angles) 
    ys = center[2] .+ radius .* sin.(angles) 
    return xs, ys
end
function click_plot_check(state::Balls_State, position)
    if((position[1]^2+position[2]^2)>0.81)
        return
    end
    if state.control % 2 == 1
        position1 = Observable(Vector{Float32}(position))
        push!(center_vec, position1)
        new_circle = @lift generate_circle($position1, 0.1, 30)
        push!(circle_points, new_circle)
        push!(xs, @lift $new_circle[1])
        push!(ys, @lift $new_circle[2])
        poly!(ax, xs[end], ys[end]; color=:black) 
        state.position_of_center = position  
        state.control += 1                 
    else
        state.alfa = atan(position[2] - state.position_of_center[2] , position[1] - state.position_of_center[1])
    end
end
function button_click_check(state::Balls_State, balls, mass_box, force_box)
    if !isnothing(mass_box.stored_string[]) && !isnothing(force_box.stored_string[])
        mass = tryparse(Float32, mass_box.stored_string[])
        force = tryparse(Float32, force_box.stored_string[])
        vector_force = [force * cos(state.alfa), force * sin(state.alfa)]       
        new_ball = Ball(mass, state.position_of_center, vector_force)
        push!(balls, new_ball)
        push!(state_vec,state.position_of_center[1])
        push!(state_vec,state.position_of_center[2])
        push!(state_vec,vector_force[1])
        push!(state_vec,vector_force[2])
        state.control += 1
        mass_box.stored_string= nothing
        force_box.stored_string= nothing
    end
end
function DiffEq(state,input,new_ball::Ball)
    Matrix_A[1, 3] = 1.0f0
    Matrix_A[2, 4] = 1.0f0
    Matrix_A[3, 3] = -c[] / new_ball.mass
    Matrix_A[4, 4] = -c[] / new_ball.mass 
    Matrix_B[3,1]=1.0f0/new_ball.mass
    Matrix_B[4,2]=1.0f0/new_ball.mass

    k1 = (Matrix_A*state + Matrix_B*input)*step
    k2 = (Matrix_A*(state + 0.5*k1) +  Matrix_B*input)*step
    k3 = (Matrix_A*(state + 0.5*k2) +  Matrix_B*input)*step
    k4 = (Matrix_A*(state +  k3) + Matrix_B*input)*step
    new_state=state + (1 / 6.0f0)*(k1 + 2 * k2 + 2 * k3 + k4)

    return new_state
end
function check_collison(ball::Ball)
    if(sqrt(ball.position[1]^2+ball.position[2]^2)>=0.9)
        return true
    else 
        return false
    end
end
center_vec=Vector{Observable{Vector{Float32}}}()
circle_points = Vector{Observable{Tuple{Vector{Float32}, Vector{Float32}}}}()
xs = Vector{Observable{Vector{Float32}}}()
ys = Vector{Observable{Vector{Float32}}}()
create_state = init_state()
fig = Figure(size=(1920, 1080))
display(fig)
ax = Axis(fig[1, 1],aspect=DataAspect())
Makie.deactivate_interaction!(ax,:rectanglezoom)
c=Observable(0.0f0)
e=Observable(0.0f0)
lines!(ax,Circle(Point2f0(0.0,0.0), 1.0), color = :black, linewidth = 10)
buttons_grid= fig[1,2]= GridLayout(tellwidth=false,tellheight=false)
mass_box = Textbox(buttons_grid[1,1],placeholder="Mass",reset_on_defocus=true, height=30,width=150)
force_box = Textbox(buttons_grid[2,1],placeholder="Force",reset_on_defocus=true,height=30,width=150)
submit = Button(buttons_grid[3,1], label="Submit",height=30,width=150)
label_air=Label(buttons_grid[4, 1], "c= $(c[])")
slider_air=Slider(buttons_grid[5,1], range=0:0.01:0.99,startvalue=0.0)
label_e=Label(buttons_grid[6, 1], "e= $(e[])")
slider_e=Slider(buttons_grid[7,1], range=0:0.01:0.99,startvalue=0.0)
Label(buttons_grid[8, 1], "Run Simulation")
toggle=Toggle(buttons_grid[9,1],height=30,width=50)

spoint = select_point(ax.scene)
balls=Vector{Ball}()
on(spoint) do position 
    if !toggle.active[]
        click_plot_check(create_state,position)
    end
end
on(slider_air.value) do new_c
    c[]=new_c
    label_air.text= "c = $(c[])"
end
on(slider_e.value) do new_e
    e[]=new_e
    label_e.text= "e = $(e[])"
end
on(submit.clicks) do g
    if !toggle.active[]
        button_click_check(create_state, balls, mass_box, force_box)
    end
end


Matrix_A = zeros(Float32, 4, 4)  # 4x4 matrix of zeros
Matrix_B = zeros(Float32, 4, 2)  # 4x2 matrix of zeros
collision_bool=false
global g=9.81f0
global step=0.005f0
state_vec=Vector{Float32}()
new_state_vec=Vector{Float32}()
input=Vector{Float32}([0.0f0,0.0f0])
# @async begin 
#     while true
#         if !isopen(fig.scene) 
#             break
#         elseif toggle.active[]
#             println("Jestem")
#             input=[balls[1].force[1],balls[1].force[2]]
#             new_state_vec=DiffEq(state_vec,input,balls[1])
#             println("Jestem")
#             state_vec=new_state_vec
#             center_vec[1][] = [state_vec[1],state_vec[2]]  # Increment the position
#             sleep(0.001)  # Controls animation speed
#         end
#         sleep(0.01)
#     end
# end
for i in 1:1000
    if i==1
        input=[balls[1].force[1]/10,balls[1].force[2]/10-g*balls[1].mass]     
    else
        input=[0.0f0,-g*balls[1].mass]
    end   
    new_state_vec=DiffEq(state_vec,input,balls[1])
    state_vec=new_state_vec
    if i%2==0
        center_vec[1][] = [state_vec[1],state_vec[2]]
    end
    balls[1].position=[state_vec[1],state_vec[2]]
    println("Speed  Vx=$(state_vec[3])  Vy=$(state_vec[4])")
    if check_collison(balls[1]) && collision_bool
        continue
    elseif !check_collison(balls[1]) && collision_bool
        collision_bool=false
    elseif check_collison(balls[1])
        println("Zderzenie!!!!!")
        dist = sqrt(balls[1].position[1]^2+balls[1].position[2]^2)
        normal_x = balls[1].position[1] / dist
        normal_y = balls[1].position[2] / dist
        v_normal = state_vec[3] * normal_x + state_vec[4] * normal_y
        if v_normal > 0
            state_vec[3] =  (state_vec[3]-(2) * v_normal * normal_x)*e[] 
            state_vec[4] =  (state_vec[4]-(2) * v_normal * normal_y)*e[]
        end
        #state_vec[4]=-state_vec[4]*e[]
        collision_bool=true
    end
    sleep(step)
end
# on(toggle.active) do y
#     @async while toggle.active[]
#         if !isopen(fig.scene) 
#             break
#         end
#         println("Jestem")
#         input=[balls[1].force[1],balls[1].force[2]]
#         new_state_vec=DiffEq(state_vec,input,balls[1])
#         state_vec=new_state_vec
#         center_vec[1][] = [state_vec[1],state_vec[2]]  # Increment the position
#         sleep(0.001)  # Controls animation speed
#     end
# end

# for _ in 1:Threads.nthreads()
#     @spawn begin
#         while true
#             if !isopen(fig.scene)
#                 break
#             elseif toggle.active[]
#                 input = balls[1].force
#                 new_state_vec = DiffEq(state_vec,input, balls[1])
#                 state_vec = new_state_vec
#                 center_vec[1][] = [state_vec[1], state_vec[2]]  # Increment the position
#                 sleep(0.001)  # Controls animation speed
#             end
#             sleep(0.01)
#         end
#     end
# end
