using GLMakie
using Observables
using GeometryBasics
using LinearAlgebra

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
    return xs,ys
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
        push!(no_bounces,false)
        push!(during_coll,false)
        push!(initial_energy,mass*g*(ball_state.position_of_center[2]+1.0)+mass*((vector_force[1]/mass)^2+(vector_force[2]/mass)^2)/2)
        push!(counter_wall,0)
        push!(coll_time,0)
        ball_state.control += 1
        mass_box.stored_string= nothing
        force_box.stored_string= nothing
    end
end
function DiffEq(state,input,new_ball::Ball,step::Float32)
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
    if(sqrt(ball.position[1]^2+ball.position[2]^2)>0.9)
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
function calculate_collision_with_ball(ball1::Ball, ball2::Ball,state1::Vector{Float32},state2::Vector{Float32})
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
function collision_with_wall(ball1::Ball,state1::Vector{Float32})
    dist = sqrt(ball1.position[1]^2+ball1.position[2]^2)
    normal_x = ball1.position[1] / dist
    normal_y = ball1.position[2] / dist 
    v_normal = state1[3] * normal_x + state1[4] * normal_y
    if v_normal > 0
        state1[3] =  (state1[3]-(2) * v_normal * normal_x)*e[] 
        state1[4] =  (state1[4]-(2) * v_normal * normal_y)*e[]
    end 
end
function deformation_with_wall(balls::Ball,kinetic_energy::Float32,counter::Int,center::Observable,new_circle::Observable,xs::Observable,ys::Observable,initial_energy::Float32)  
    alfa=kinetic_energy/initial_energy
    deformation_per=(1-e[])*alfa/4
    dist = sqrt(balls.position[1]^2 + balls.position[2]^2)
    if dist > 0.9
        # Wektor normalny
        normal_x = balls.position[1] / dist
        normal_y = balls.position[2] / dist
    
        if counter==9
            balls.position[1] -= (dist - 0.9) * normal_x
            balls.position[2] -= (dist - 0.9) * normal_y 
        elseif counter>4
            balls.position[1] -= 0.2*deformation_per * normal_x
            balls.position[2] -= 0.2*deformation_per * normal_y  
        else
            balls.position[1] += 0.2*deformation_per * normal_x
            balls.position[2] += 0.2*deformation_per * normal_y     
        end 
        center[]=[balls.position[1],balls.position[2]]
    end
    points = new_circle[]
    points_x=points[1]
    points_y=points[2]
    combined = [[points_x[i], points_y[i]] for i in 1:30] 
    temp_x=Vector{Float32}()
    temp_y=Vector{Float32}()
    for point in combined       
        dist = sqrt(point[1]^2 + point[2]^2)    
        if dist>1.0
            normal_x = point[1] / dist
            normal_y = point[2] / dist
            point[1]-= (dist - 1) * normal_x
            point[2]-= (dist - 1) * normal_y     
        end
        push!(temp_x,point[1])
        push!(temp_y,point[2])
    end

    xs[]= temp_x
    ys[]= temp_y
    
end
# function deformation_with_ball(ball1::Ball,ball2::Ball,state1::Vector{Float32},state2::Vector{Float32},new_circle1::Observable,new_circle2::Observable,xs1::Observable,xs2::Observable,ys1::Observable,ys2::Observable,center1::Observable,center2::Observable)
#     normal = normalize(ball2.position - ball1.position)
#     #dist = norm(ball2.position - ball1.position)
#     #target_dist = 0.08  # 40% średnicy
#     offset = 0.06/2
    
#     ball1.position +=  offset * normal
#     ball2.position -=  offset * normal
#     center1[]=ball1.position
#     center2[]=ball2.position
#     println(norm(ball2.position - ball1.position))
#     points = new_circle1[]
#     points_x=points[1]
#     points_y=points[2]
#     combined = [[points_x[i], points_y[i]] for i in 1:30] 
#     temp_x=Vector{Float32}()
#     temp_y=Vector{Float32}()
#     for point in combined       
       
#         push!(temp_x,point[1])
#         push!(temp_y,point[2])
#     end

#     points2 = new_circle2[]
#     points_x2=points2[1]
#     points_y2=points2[2]
#     combined2 = [[points_x2[i], points_y2[i]] for i in 1:30] 
#     temp_x2=Vector{Float32}()
#     temp_y2=Vector{Float32}()
#     for point in combined2       
       
#         push!(temp_x2,point[1])
#         push!(temp_y2,point[2])
#     end
    
#     xs1[]= temp_x
#     ys1[]= temp_y
#     xs2[]= temp_x2
#     ys2[]= temp_y2
    
# end
function correction(state1::Vector{Float32},state2::Vector{Float32})
    dist = sqrt((state1[1]-state2[1])^2 + (state1[2]-state2[2])^2)
    normal_x = (state2[1] - state1[1])/ dist
    normal_y = (state2[2] - state1[2]) / dist

    state1[1] -= (0.201-dist) * normal_x/2
    state1[2] -= (0.201-dist) * normal_y/2
    state2[1] += (0.201-dist) * normal_x/2
    state2[2] += (0.201-dist) * normal_y/2
end
function start() 
    for i in 1:1000
        if !isopen(fig.scene) || !isrunning[]
            return
        end       
        for j in 1:length(balls)            
            if no_bounces[j]
                input_vec[j].=[0.0f0,0.0f0]
                state_vec[j][3],state_vec[j][4]=[0.0f0,0.0f0] 
            else
                input_vec[j].=[0.0f0,-g*balls[j].mass]            
            end  
 
            if check_collison(balls[j]) && !no_bounces[j]        
                if e[]!=1.0f0
                    coll_time[j]+=1 
                    kinetic_energy=balls[j].mass*(state_vec[j][3]^2+state_vec[j][4]^2)/2  
                    if abs(state_vec[j][3])>0.9 || abs(state_vec[j][4])>0.9  deformation_with_wall(balls[j],kinetic_energy,coll_time[j],center_vec[j],circle_points[j],xs[j],ys[j],initial_energy[j]) end     
                end               
                if coll_time[j]<9 && e[]!=1.0f0 
                    continue
                else
                    coll_time[j]=0
                    if no_bounces[j]
                         continue
                    end
                    
                    if abs(state_vec[j][3]) <0.1 && abs(state_vec[j][4]) <0.1 && !no_bounces[j] && state_vec[j][2]<0
                        no_bounces[j]=true
                        continue
                    end  
                    collision_with_wall(balls[j],state_vec[j])
                    if e[]!=1.0f0 counter_wall[j]+=1 end
                end
            end
            for coll in j+1:length(balls)
                #dist=sqrt((balls[j].position[1]-balls[coll].position[1])^2+(balls[j].position[2]-balls[coll].position[2])^2)
                collision_occ=check_collision_balls(j,coll)
                         
                if collision_occ
                    correction(state_vec[j],state_vec[coll])                 
                end
                # if collision_occ && e[]!=1
                #     deformation_with_ball(balls[j],balls[coll],state_vec[j],state_vec[coll],circle_points[j],circle_points[coll],xs[j],xs[coll],ys[j],ys[coll],center_vec[j],center_vec[coll])
                #     poly!(ax, xs[j], ys[j].+0.5; color=:black) 
                #     poly!(ax, xs[coll], ys[coll].+0.5; color=:black) 
                #     sleep(10)
                # else
                if collision_occ && no_bounces[j] && no_bounces[coll] 
                    continue
                elseif collision_occ && (no_bounces[j] || no_bounces[coll] ) 
                    if no_bounces[j]
                        normal = normalize(balls[coll].position - balls[j].position)
                        v1=[state_vec[coll][3],state_vec[coll][4]]
                        v1n = dot(v1, normal)
                        v_new = v1 .- (1 + e[]) * v1n * normal
                        state_vec[coll][3],state_vec[coll][4]=v_new
                    elseif no_bounces[coll]
                        normal = normalize(balls[j].position - balls[coll].position)
                        v1=[state_vec[j][3],state_vec[j][4]]
                        v1n = dot(v1, normal)
                        v_new = v1 .- (1 + e[]) * v1n * normal
                        state_vec[j][3],state_vec[j][4]=v_new
                    end

                elseif collision_occ 
                    calculate_collision_with_ball(balls[j],balls[coll],state_vec[j],state_vec[coll])       
                end
                #end
            end
            
        
            state_vec[j].=DiffEq(state_vec[j],input_vec[j],balls[j],step)

            if sqrt(state_vec[j][1]^2+state_vec[j][2]^2)>0.9 
                normal_x = state_vec[j][1] / sqrt(state_vec[j][1]^2+state_vec[j][2]^2)
                normal_y = state_vec[j][2] / sqrt(state_vec[j][1]^2+state_vec[j][2]^2)
                state_vec[j][1]-= (sqrt(state_vec[j][1]^2+state_vec[j][2]^2) - 0.9) * normal_x
                state_vec[j][2]-= (sqrt(state_vec[j][1]^2+state_vec[j][2]^2) - 0.9) * normal_y    
            end
            center_vec[j][] = [state_vec[j][1],state_vec[j][2]]
            balls[j].position=[state_vec[j][1],state_vec[j][2]]  
            
        end 
        sleep(step)
    end  
           
end

const g=9.81f0
const step=0.002f0
Matrix_A = zeros(Float32, 4, 4)  # 4x4 matrix of zeros
Matrix_B = zeros(Float32, 4, 2)  # 4x2 matrix of zeros
isrunning=Observable(false)
c=Observable(0.0f0)
e=Observable(0.0f0)

balls=Vector{Ball}()
state_vec=Vector{Vector{Float32}}()
input_vec=Vector{Vector{Float32}}()
no_bounces=Vector{Bool}()
during_coll=Vector{Bool}()
counter_wall=Vector{Int}()
initial_energy=Vector{Float32}()
coll_time=Vector{Int}()
center_vec=Vector{Observable{Vector{Float32}}}()
circle_points = Vector{Observable{Tuple{Vector{Float32}, Vector{Float32}}}}()
xs = Vector{Observable{Vector{Float32}}}()
ys = Vector{Observable{Vector{Float32}}}()

create_state = init_state()
fig = Figure(size=(1920, 1080))
display(fig)
ax = Axis(fig[1, 1],aspect=DataAspect())
Makie.deactivate_interaction!(ax,:rectanglezoom)
lines!(ax,Circle(Point2f0(0.0,0.0), 1.0), color = :black, linewidth = 10)
buttons_grid= fig[1,2]= GridLayout(tellwidth=false,tellheight=false)
mass_box = Textbox(buttons_grid[1,1],placeholder="Mass",reset_on_defocus=true, height=30,width=150)
force_box = Textbox(buttons_grid[2,1],placeholder="Force",reset_on_defocus=true,height=30,width=150)
submit = Button(buttons_grid[3,1], label="Submit",height=30,width=150)
label_air=Label(buttons_grid[4, 1], "c= $(c[])")
slider_air=Slider(buttons_grid[5,1], range=0:0.01:5,startvalue=0.0)
label_e=Label(buttons_grid[6, 1], "e= $(e[])")
slider_e=Slider(buttons_grid[7,1], range=0:0.01:1,startvalue=0.0)
toggle=Button(buttons_grid[8,1],label="Run Simulation",height=30,width=150)
reset_button=Button(buttons_grid[9,1],label="Reset",height=30,width=150)
spoint = select_point(ax.scene)


on(reset_button.clicks) do reset
    isrunning[]=false
    empty!(balls)
    empty!(state_vec)
    empty!(input_vec)
    empty!(no_bounces)
    empty!(during_coll)
    empty!(counter_wall)
    empty!(initial_energy)
    empty!(coll_time)
    empty!(center_vec)
    empty!(circle_points)
    empty!(xs)
    empty!(ys)
    empty!(ax)
    lines!(ax,Circle(Point2f0(0.0,0.0), 1.0), color = :black, linewidth = 10)
end
on(spoint) do position 
    if !isrunning[]
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
    if !isrunning[]
        button_click_check(create_state, balls, mass_box, force_box)
    end
end
on(toggle.clicks) do click 
    if !isrunning[]
        isrunning[] = true
        toggle.label = "Running..."
        @async begin
            if !isopen(fig.scene)
                return
            end
            start()
            isrunning[] = false
            toggle.label = "Run"
        end
    end  
end