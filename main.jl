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
# mutable struct state
#     p_x::Float32
#     p_y::Float32
#     v_x::Float32
#     v_y::Float32
# end
# mutable struct input
#     f_x::Float32
#     f_y::Float32
# end
function init_state()
    init=Balls_State(1,[0.0,0.0],0.0)
end
function generate_circle(center::Vector{Float32}, radius::Float64, n_points::Int)
    angles = LinRange(0, 2pi, n_points) 
    xs = center[1] .+ radius .* cos.(angles) 
    ys = center[2] .+ radius .* sin.(angles) 
    return xs, ys
end
function click_plot_check(ball_state::Balls_State, position)
    if((position[1]^2+position[2]^2)>0.8)
        return
    end
    if ball_state.control % 2 == 1
        position1 = Observable(Vector{Float32}(position))
        push!(center_vec, position1)
        new_circle = @lift generate_circle($position1, 0.1, 30)
        push!(circle_points, new_circle)
        push!(xs, @lift $new_circle[1])
        push!(ys, @lift $new_circle[2])
        poly!(ax, xs[end], ys[end]; color=:black) 
        ball_state.position_of_center = position  
        ball_state.control += 1                 
    else
        ball_state.alfa = atan(position[2] - ball_state.position_of_center[2] , position[1] - ball_state.position_of_center[1])
    end
end
function button_click_check(ball_state::Balls_State, balls, mass_box, force_box)
    if !isnothing(mass_box.stored_string[]) && !isnothing(force_box.stored_string[])
        mass = tryparse(Float32, mass_box.stored_string[])
        force = tryparse(Float32, force_box.stored_string[])
        vector_force = [force * cos(ball_state.alfa), force * sin(ball_state.alfa)]       
        new_ball = Ball(mass, ball_state.position_of_center, vector_force)
        push!(balls, new_ball)
        push!(state_vec,[ball_state.position_of_center[1],ball_state.position_of_center[2],vector_force[1]/mass,vector_force[2]/mass])
        push!(input_vec,[0.0f0,0.0f0])
        ball_state.control += 1
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
function check_collision_balls(i,j)
    if(sqrt((balls[i].position[1]-balls[j].position[1])^2+(balls[i].position[2]-balls[j].position[2])^2)<=0.2)
        return true
    else 
        return false
    end  
end
function calculate_collision(ball1::Ball, ball2::Ball,state1::Vector{Float32},state2::Vector{Float32})
 # Wektor normalny między środkami
    normal = normalize(ball2.position - ball1.position)
    tangent = [-normal[2], normal[1]]  # Wektor styczny dla osi x i y
    v1=[state1[3],state1[4]]
    v2=[state2[3],state2[4]]
    # Prędkości normalne i styczne przed zderzeniem (skalarne)
    v1n = dot(v1, normal)
    v2n = dot(v2, normal)
    v1t = dot(v1, tangent)
    v2t = dot(v2, tangent)

    # Prędkości normalne po zderzeniu (z uwzględnieniem niesprężystości)
    v1n_after = ((ball1.mass - e[] * ball2.mass) * v1n + (1 + e[]) * ball2.mass * v2n) / (ball1.mass + ball2.mass)
    v2n_after = ((ball2.mass - e[] * ball1.mass) * v2n + (1 + e[]) * ball1.mass * v1n) / (ball1.mass + ball2.mass)

    # Prędkości styczne pozostają niezmienione (brak tarcia między kulkami)
    v1t_after = v1t
    v2t_after = v2t

    # Nowe prędkości dla osi x i y po zderzeniu
    state1[3] = v1n_after * normal[1] + v1t_after * tangent[1]
    state1[4] = v1n_after * normal[2] + v1t_after * tangent[2]
    state2[3] = v2n_after * normal[1] + v2t_after * tangent[1]
    state2[4] = v2n_after * normal[2] + v2t_after * tangent[2]
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
slider_air=Slider(buttons_grid[5,1], range=0:0.01:1,startvalue=0.0)
label_e=Label(buttons_grid[6, 1], "e= $(e[])")
slider_e=Slider(buttons_grid[7,1], range=0:0.01:1,startvalue=0.0)
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
on(submit.clicks) do check
    if !toggle.active[]
        button_click_check(create_state, balls, mass_box, force_box)
    end
end

global no_bounces=false
Matrix_A = zeros(Float32, 4, 4)  # 4x4 matrix of zeros
Matrix_B = zeros(Float32, 4, 2)  # 4x2 matrix of zeros
const friction_coefficient=0.2f0
const g=9.81f0
const step=0.005f0
state_vec=Vector{Vector{Float32}}()
input_vec=Vector{Vector{Float32}}()
if toggle.active[]
    for i in 1:1000
        for j in 1:length(balls)
            if i==1
                input_vec[j].=[balls[j].force[1],balls[j].force[2]-g*balls[j].mass]
            # elseif no_bounces
            #     force_normal=balls[1].mass*g/cos(atan(state_vec[2],state_vec[1]))
            #     input.=[g*balls[1].mass*sin(atan(state_vec[2]/state_vec[1])),-g*balls[1].mass]
            #     #input.=[force_normal*cos(atan(state_vec[2],state_vec[1]))*(1-friction_coefficient),-g*balls[1].mass+force_normal*sin(atan(state_vec[2],state_vec[1]))*(1-friction_coefficient)]
            else
                input_vec[j].=[0.0f0,-g*balls[j].mass]
            end  

            state_vec[j].=DiffEq(state_vec[j],input_vec[j],balls[j])

            if i%2==0 && check_collison(balls[j]) && no_bounces 
                state_vec[j][2]-=(state_vec[j][2]+sqrt(0.81-state_vec[j][1]^2))
                center_vec[j][] = [state_vec[j][1],state_vec[j][2]]
            elseif i%2==0
                if state_vec[j][1]^2+state_vec[j][2]^2 >= 1
                    state_vec[j][1]-=(state_vec[j][1]+sqrt(0.81-state_vec[j][2]^2))
                    state_vec[j][2]-=(state_vec[j][2]+sqrt(0.81-state_vec[j][1]^2))
                else
                    center_vec[j][] = [state_vec[j][1],state_vec[j][2]] 
                end    
            end  
            
            balls[j].position=[state_vec[j][1],state_vec[j][2]]
            #println("Speed  Vx=$(state_vec[3])  Vy=$(state_vec[4])")   
            if check_collison(balls[j])
                if no_bounces
                 continue
                end
                if abs(state_vec[j][3]) <0.1 && abs(state_vec[j][4]) <0.1 && !no_bounces && state_vec[j][2]<0
                    println("no bounces")
                    global no_bounces=true
                    continue
                end  
                println("Zderzenie!!!!!")
                dist = sqrt(balls[j].position[1]^2+balls[j].position[2]^2)
                normal_x = balls[j].position[1] / dist
                normal_y = balls[j].position[2] / dist 
                v_normal = state_vec[j][3] * normal_x + state_vec[j][4] * normal_y
                if v_normal > 0
                    state_vec[j][3] =  (state_vec[j][3]-(2) * v_normal * normal_x)*e[] 
                    state_vec[j][4] =  (state_vec[j][4]-(2) * v_normal * normal_y)*e[]
                end   
            end
            for coll in j+1:length(balls)
                if check_collision_balls(j,coll)
                    println("Jestem")
                    calculate_collision(balls[j],balls[coll],state_vec[j],state_vec[coll])
                end
            end
            
        end
        sleep(step)
    end
end
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
