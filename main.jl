using GLMakie
using Observables
using GeometryBasics
using LinearAlgebra

mutable struct Ball
    mass::Float32
    position::Vector{Float32}
    force::Vector{Float32}
    velocity::Vector{Float32}
end
mutable struct Balls_State
    control::Int
    position_of_center::Vector{Float32}
    alfa::Float32
end
function generate_circle(center::Vector{Float32}, radius::Float64, n_points::Int)
    angles = LinRange(0, 2pi, n_points) 
    x_circle = center[1] .+ radius .* cos.(angles) 
    y_circle = center[2] .+ radius .* sin.(angles) 
    return x_circle,y_circle
end
function click_plot_check(ball_state::Balls_State, position)
    if((position[1]^2+position[2]^2)>0.8) return end
  
    if ball_state.control % 2 == 1
        position1 = Observable(Vector{Float32}(position))
        push!(center_vec, position1)
        new_circle = @lift generate_circle($position1, 0.1, 30)
        push!(circle_points, new_circle)
        push!(x_circle, @lift $new_circle[1])
        push!(y_circle, @lift $new_circle[2])
        poly!(ax, x_circle[end], y_circle[end]; color=:black) 
        ball_state.position_of_center = position  
        ball_state.control += 1                
    else
        ball_state.alfa = atan(position[2] - ball_state.position_of_center[2] , position[1] - ball_state.position_of_center[1])
    end
end
function button_click_check(ball_state::Balls_State, balls::Vector{Ball}, mass_box, force_box)
    if isnothing(mass_box.stored_string[]) || isnothing(force_box.stored_string[]) || mass_box.stored_string[]=="0" return end        
    mass = tryparse(Float32, mass_box.stored_string[])
    force = tryparse(Float32, force_box.stored_string[])
    if isnothing(mass) || isnothing(force) return end
    
    vector_force = [min(force,mass*max_velocity) * cos(ball_state.alfa), min(force,mass*max_velocity) * sin(ball_state.alfa)]       
    new_ball = Ball(mass, ball_state.position_of_center, vector_force,[vector_force[1]/mass,vector_force[2]/mass])
    push!(balls, new_ball)
    push!(state_vec,[ball_state.position_of_center[1],ball_state.position_of_center[2],vector_force[1]/mass,vector_force[2]/mass])
    push!(input_vec,[0.0f0,0.0f0])
    push!(no_bounces,false)
    push!(initial_velocity,max_velocity)
    push!(coll_time,0)
    push!(coll_time_ball,0)
    ball_state.control += 1
    mass_box.stored_string= nothing
    force_box.stored_string= nothing
end
function DiffEq(state,input,new_ball,step)
    Matrix_A[1, 3] = 1.0f0
    Matrix_A[2, 4] = 1.0f0
    Matrix_A[3, 3] = -c[] / new_ball.mass
    Matrix_A[4, 4] = -c[]/ new_ball.mass 
    Matrix_B[3,1]=1.0f0/new_ball.mass
    Matrix_B[4,2]=1.0f0/new_ball.mass

    k1 = (Matrix_A*state + Matrix_B*input)*step
    k2 = (Matrix_A*(state + 0.5*k1) +  Matrix_B*input)*step
    k3 = (Matrix_A*(state + 0.5*k2) +  Matrix_B*input)*step
    k4 = (Matrix_A*(state +  k3) + Matrix_B*input)*step
    
    new_state=state + (1 / 6.0f0)*(k1 + 2 * k2 + 2 * k3 + k4)
    return new_state
end
function check_collison(j)
    return (sqrt(balls[j].position[1]^2+balls[j].position[2]^2)>0.91) ? true : false 
end
function check_collision_balls(i,j)
   return (sqrt((balls[i].position[1]-balls[j].position[1])^2+(balls[i].position[2]-balls[j].position[2])^2)<=0.2) ? true : false 
end
function collision_with_ball(ball1::Ball, ball2::Ball,state1::Vector{Float32},state2::Vector{Float32})
    normal = normalize(ball2.position - ball1.position)
    tangent = [-normal[2], normal[1]] 

    v1n = dot([state1[3],state1[4]], normal)
    v2n = dot([state2[3],state2[4]], normal)
    v1t = dot([state1[3],state1[4]], tangent)
    v2t = dot([state2[3],state2[4]], tangent)

    v1n_after = ((ball1.mass - e[] * ball2.mass) * v1n + (1 + e[]) * ball2.mass * v2n) / (ball1.mass + ball2.mass)
    v2n_after = ((ball2.mass - e[] * ball1.mass) * v2n + (1 + e[]) * ball1.mass * v1n) / (ball1.mass + ball2.mass)

    state1[3] = v1n_after * normal[1] + v1t * tangent[1]
    state1[4] = v1n_after * normal[2] + v1t * tangent[2]
    state2[3] = v2n_after * normal[1] + v2t * tangent[1]
    state2[4] = v2n_after * normal[2] + v2t * tangent[2]
end
function collision_with_no_bouncing_ball(ball1::Ball,ball2::Ball,state::Vector{Float32})
    normal = normalize(ball2.position - ball1.position)
    v1=[state[3],state[4]]
    v1n = dot(v1, normal)
    v_new = v1 .- (1 + e[]) * v1n * normal 
    state[3],state[4]=v_new
end
function collision_with_wall(ball1::Ball,state1::Vector{Float32}) 
    normal=normalize(ball1.position)
    v_normal = state1[3] * normal[1] + state1[4] * normal[2]
    if v_normal > 0 state1[3],state1[4] = [ (state1[3]-(2) * v_normal * normal[1])*e[],(state1[4]-(2) * v_normal * normal[2])*e[]] end
end
function deformation_with_wall(balls::Ball,counter::Int,center::Observable,new_circle::Observable,x_circle::Observable,y_circle::Observable,initial_velocity::Float32)  
    alfa=min(norm(balls.velocity)/initial_velocity,2.0)
    deformation_per=k[]*alfa
    if deformation_per<0.02 return deformation_per end

    dist = sqrt(balls.position[1]^2 + balls.position[2]^2)
    if dist > 0.9
        normal_x = balls.position[1] / dist
        normal_y = balls.position[2] / dist
        if counter==6
            balls.position[1] -= (dist - 0.9) * normal_x
            balls.position[2] -= (dist - 0.9) * normal_y 
        elseif counter>3
            balls.position[1] -= 0.2*deformation_per * normal_x/3
            balls.position[2] -= 0.2*deformation_per * normal_y/3
        else
            balls.position[1] += 0.2*deformation_per * normal_x/3
            balls.position[2] += 0.2*deformation_per * normal_y/3    
        end 
        center[]=[balls.position[1],balls.position[2]]
    end
    points = new_circle[]
    combined_points = [[points[1][i], points[2][i]] for i in 1:30] 
    temp_x=Vector{Float32}()
    temp_y=Vector{Float32}()
    for point in combined_points       
        dist = sqrt(point[1]^2 + point[2]^2)    
        if dist>1.0
            normal=normalize([point[1],point[2]])  
            point -= (dist-1)*normal  
        end
        push!(temp_x,point[1])
        push!(temp_y,point[2])
    end
    x_circle[]= temp_x
    y_circle[]= temp_y

    return 1.0
end
function deformation_with_ball(ball1::Ball,ball2::Ball,time::Int,initial_velocity1,initial_velocity2,center1,center2)
    alfa1=min(norm(ball1.velocity)/initial_velocity1,2.0)
    alfa2=min(norm(ball2.velocity)/initial_velocity2,2.0)
    deformation_per1=(ball2.mass/(ball2.mass+ball1.mass))*k[]*alfa1
    deformation_per2=(ball1.mass/(ball2.mass+ball1.mass))*k[]*alfa2
    if deformation_per1+deformation_per2 <0.02
        return deformation_per1+deformation_per2
    end
    
    normal = normalize(ball2.position - ball1.position)
    if time < 4
        ball1.position += 0.2*deformation_per1 * normal/3
        ball2.position -= 0.2*deformation_per2 * normal/3
    else
        ball1.position -= 0.2*deformation_per1* normal/3
        ball2.position += 0.2*deformation_per2 * normal/3
    end
    center1[]=[ball1.position[1],ball1.position[2]]
    center2[]=[ball2.position[1],ball2.position[2]]
    
    return 1.0
end
function correction(state1::Vector{Float32},state2::Vector{Float32})
    dist = sqrt((state1[1]-state2[1])^2 + (state1[2]-state2[2])^2)
    normal=normalize([state2[1],state2[2]] - [state1[1],state1[2]])
    state1[1] -= (0.201-dist) * normal[1]/2
    state1[2] -= (0.201-dist) * normal[2]/2
    state2[1] += (0.201-dist) * normal[1]/2
    state2[2] += (0.201-dist) * normal[2]/2
end
function start() 
    while true
        if !isopen(fig.scene) || !isrunning[] return end       
        for j in 1:length(balls)            
            if no_bounces[j]
                input_vec[j].=[0.0f0,0.0f0]
                state_vec[j][3],state_vec[j][4]=[0.0f0,0.0f0] 
            else
                input_vec[j].=[0.0f0,-g*balls[j].mass]            
            end  
 
            if check_collison(j) && !no_bounces[j]           
                no_deformation=false
                coll_time[j]+=1  

                if deformation_with_wall(balls[j],coll_time[j],center_vec[j],circle_points[j],x_circle[j],y_circle[j],initial_velocity[j]) < 0.02 no_deformation=true end 
                
                if coll_time[j]<6  && !no_deformation
                    continue
                elseif abs(state_vec[j][3]) <0.1 && abs(state_vec[j][4]) <0.3  && !no_bounces[j] && state_vec[j][2]<0
                    no_bounces[j]=true
                    balls[j].mass=9999999.0f0
                else
                    coll_time[j]=0
                    collision_with_wall(balls[j],state_vec[j])
                end  
            end
            for coll in j+1:length(balls)
                collision_occ=check_collision_balls(j,coll)
               
                if collision_occ 
                    no_deformation_ball=false
                    coll_time_ball[j]+=1 
                    coll_time_ball[coll]+=1 
                    if no_bounces[j] && no_bounces[coll] correction(state_vec[j],state_vec[coll]) end

                    if deformation_with_ball(balls[j],balls[coll],coll_time_ball[j],initial_velocity[j],initial_velocity[coll],center_vec[j],center_vec[coll]) < 0.02 no_deformation_ball=true end
                    
                    if  coll_time_ball[j]<6 && coll_time_ball[coll]<6 && !no_deformation_ball continue end

                    coll_time_ball[j]=0
                    coll_time_ball[coll]=0
                    correction(state_vec[j],state_vec[coll])
                    if no_bounces[j] || no_bounces[coll]  
                        if no_bounces[j]
                            collision_with_no_bouncing_ball(balls[j],balls[coll],state_vec[coll])
                        elseif no_bounces[coll]
                            collision_with_no_bouncing_ball(balls[coll],balls[j],state_vec[j])
                        end
                    else
                        collision_with_ball(balls[j],balls[coll],state_vec[j],state_vec[coll])       
                    end
                end
            end
           
            if coll_time_ball[j]!=0
                continue
            else
                state_vec[j].=DiffEq(state_vec[j],input_vec[j],balls[j],step)  
            end
            
            if sqrt(state_vec[j][1]^2+state_vec[j][2]^2)> 0.91
                dist=sqrt(state_vec[j][1]^2+state_vec[j][2]^2)
                normal=normalize([state_vec[j][1],state_vec[j][2]])
                state_vec[j][1] -= (dist - 0.9101) *  normal[1]
                state_vec[j][2] -= (dist - 0.9101) *  normal[2]
            end
            center_vec[j][] = [state_vec[j][1],state_vec[j][2]]
            balls[j].position=[state_vec[j][1],state_vec[j][2]] 
            balls[j].velocity=[state_vec[j][3],state_vec[j][4]]
            
        end 
        sleep(step)
    end             
end


const g=9.81f0
const step=0.002f0
const max_velocity=15.0f0
Matrix_A = zeros(Float32, 4, 4)  # 4x4 matrix of zeros
Matrix_B = zeros(Float32, 4, 2)  # 4x2 matrix of zeros
isrunning=Observable(false)
c=Observable(0.0f0)
e=Observable(0.0f0)
k=Observable(0.0f0)

balls=Vector{Ball}()
state_vec=Vector{Vector{Float32}}()
input_vec=Vector{Vector{Float32}}()
no_bounces=Vector{Bool}()
initial_velocity=Vector{Float32}()
coll_time=Vector{Int}()
coll_time_ball=Vector{Int}()
center_vec=Vector{Observable{Vector{Float32}}}()
circle_points = Vector{Observable{Tuple{Vector{Float32},Vector{Float32}}}}()
x_circle = Vector{Observable{Vector{Float32}}}()
y_circle = Vector{Observable{Vector{Float32}}}()
create_state = Balls_State(1,[0.0,0.0],0.0)


fig = Figure(size=(1920, 1080))
ax = Axis(fig[1, 1],aspect=DataAspect())
spoint = select_point(ax.scene)
Makie.deactivate_interaction!(ax,:rectanglezoom)
lines!(ax,Circle(Point2f0(0.0,0.0), 1.0), color = :black, linewidth = 10)
buttons_grid= fig[1,2]= GridLayout(tellwidth=false,tellheight=false)
mass_box = Textbox(buttons_grid[1,1],placeholder="Mass",reset_on_defocus=true, height=30,width=150)
force_box = Textbox(buttons_grid[2,1],placeholder="Force",reset_on_defocus=true,height=30,width=150)
submit = Button(buttons_grid[3,1], label="Submit",height=30,width=150)
label_air=Label(buttons_grid[4, 1], "c= $(c[])")
slider_air=Slider(buttons_grid[5,1], range=0:0.01:3,startvalue=0.0)
label_e=Label(buttons_grid[6, 1], "e= $(e[])")
slider_e=Slider(buttons_grid[7,1], range=0.01:0.01:1,startvalue=0.01)
label_k=Label(buttons_grid[8, 1], "k= $(k[])")
slider_k=Slider(buttons_grid[9,1], range=0:0.01:1,startvalue=0.0)
toggle=Button(buttons_grid[10,1],label="Run Simulation",height=30,width=150)
reset_button=Button(buttons_grid[11,1],label="Reset",height=30,width=150)
display(fig)

on(reset_button.clicks) do reset
    isrunning[]=false
    empty!(balls)
    empty!(state_vec)
    empty!(input_vec)
    empty!(no_bounces)
    empty!(initial_velocity)
    empty!(coll_time)
    empty!(center_vec)
    empty!(circle_points)
    empty!(x_circle)
    empty!(y_circle)
    create_state.control=1
    create_state.position_of_center=[0.0f0,0.0f0]
    create_state.alfa=0.0f0
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
on(slider_k.value) do new_k
    k[]=new_k
    label_k.text= "k = $(k[])"
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