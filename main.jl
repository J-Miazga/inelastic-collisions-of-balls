# Informacje dotyczące fizyki zderzenie zaczerpnięte z Internetu oraz pomocy Chat GPT. Interfejs GLMakie https://docs.makie.org/v0.21/explanations/backends/glmakie
using GLMakie
using GeometryBasics
using LinearAlgebra
#struktura przechowująca dane o kulce
mutable struct Ball
    mass::Float64
    position::Vector{Float64}
    force::Vector{Float64}
    velocity::Vector{Float64}
end
#struktura przechowujaca informacje o stanie tworzenie kulek (potrzebna do wstępnego ustawienia kulek)
mutable struct Initial_State
    control::Int
    position_of_center::Vector{Float64}
    alfa::Float64
end
#deklaracja stałych, observables i wektorów
const G=9.81
const STEP=0.002
const MAX_VELOCITY=15.0
const DRAG_COEFFICIENT=0.009039275
Matrix_A = zeros(Float64, 4, 4) 
Matrix_B = zeros(Float64, 4, 2)  
is_simulating=Observable(false)
e=Observable(0.0)
k=Observable(0.0)
balls=Vector{Ball}()
state_vec=Vector{Vector{Float64}}()
input_vec=Vector{Vector{Float64}}()
no_bounces=Vector{Bool}()
initial_velocity=Vector{Float64}()
coll_time=Vector{Int}()
coll_time_ball=Vector{Int}()
center_vec=Vector{Observable{Vector{Float64}}}()
circle_points = Vector{Observable{Tuple{Vector{Float64},Vector{Float64}}}}()
x_circle = Vector{Observable{Vector{Float64}}}()
y_circle = Vector{Observable{Vector{Float64}}}()
create_state = Initial_State(1,[0.0,0.0],0.0)
#funkcja tworząca punkty będące częścią kulki z pomocą Chat GPT
function generate_circle_points(center_of_ball::Vector{Float64}, radius::Float64, number_of_points::Int)
    angle = LinRange(0, 2pi, number_of_points) 
    x_circle = center_of_ball[1] .+ radius .* cos.(angle) 
    y_circle = center_of_ball[2] .+ radius .* sin.(angle) 
    return x_circle,y_circle
end
#funkcja dodająca kulke graficznie
function click_plot_check(initial_state::Initial_State, position)
    if((position[1]^2+position[2]^2)>0.8) return end
    #uzależnienie wyświetlania kulki od observables za pomocą macro @lift
    if initial_state.control % 2 == 1
        position_ball = Observable(Vector{Float64}(position))
        push!(center_vec, position_ball)
        circle = @lift generate_circle_points($position_ball, 0.1, 30)
        push!(circle_points, circle)
        push!(x_circle, @lift $circle[1])
        push!(y_circle, @lift $circle[2])
        poly!(axis, x_circle[end], y_circle[end]; color=:black) 
        initial_state.position_of_center = position  
        initial_state.control += 1                
    else
        initial_state.alfa = atan(position[2] - initial_state.position_of_center[2] , position[1] - initial_state.position_of_center[1])
    end
end
#funkca dodająca wartości kulki do odpowiednich wektorów
function button_click_check(initial_state::Initial_State, balls::Vector{Ball}, mass_box::Textbox, force_box::Textbox)
    if isnothing(mass_box.stored_string[]) || isnothing(force_box.stored_string[]) || mass_box.stored_string[]=="0" return end        
    mass = tryparse(Float64, mass_box.stored_string[])
    force = tryparse(Float64, force_box.stored_string[])
    if isnothing(mass) || isnothing(force) return end
    
    vector_force = [min(force,mass*MAX_VELOCITY)*cos(initial_state.alfa), min(force,mass*MAX_VELOCITY)*sin(initial_state.alfa)]       
    ball = Ball(mass, initial_state.position_of_center, vector_force, [vector_force[1]/mass,vector_force[2]/mass])
    push!(balls, ball)
    push!(state_vec, [initial_state.position_of_center[1], initial_state.position_of_center[2], vector_force[1]/mass, vector_force[2]/mass])
    push!(input_vec, [0.0,0.0])
    push!(no_bounces, false)
    push!(initial_velocity, MAX_VELOCITY)
    push!(coll_time, 0)
    push!(coll_time_ball, 0)
    initial_state.control += 1
    mass_box.stored_string = nothing
    force_box.stored_string = nothing
end
#funkcja licząca równanie różniczkowe metodą Rungego-Kutty
function DiffEq(state::Vector{Float64}, input::Vector{Float64}, new_ball::Ball, STEP::Float64)
    Matrix_A[1, 3] = 1.0
    Matrix_A[2, 4] = 1.0
    Matrix_A[3, 3] = -DRAG_COEFFICIENT*sqrt(state[3]^2+state[4]^2)/new_ball.mass
    Matrix_A[4, 4] = -DRAG_COEFFICIENT*sqrt(state[3]^2+state[4]^2)/new_ball.mass 
    Matrix_B[3,1]= 1.0/new_ball.mass
    Matrix_B[4,2]= 1.0/new_ball.mass

    k1 = (Matrix_A*state + Matrix_B*input)*STEP
    k2 = (Matrix_A*(state + 0.5*k1) +  Matrix_B*input)*STEP
    k3 = (Matrix_A*(state + 0.5*k2) +  Matrix_B*input)*STEP
    k4 = (Matrix_A*(state +  k3) + Matrix_B*input)*STEP
    
    new_state=state + (1/6.0)*(k1 + 2*k2 + 2*k3 + k4)
    return new_state
end
#funkcja sprawdzająca kolizje ze ścianą
function check_collison(j::Int)
    return (norm(balls[j].position)>0.91) ? true : false 
end
#funkcja sprawdzająca kolizje ze inną kulką
function check_collision_balls(i::Int,j::Int)
   return (norm(balls[j].position - balls[i].position)<=0.2) ? true : false 
end
#funkcja licząca kolizje z inną kulką
function collision_with_ball(ball1::Ball, ball2::Ball, state1::Vector{Float64}, state2::Vector{Float64})
    normal = normalize(ball2.position - ball1.position)
    tangent = [-normal[2], normal[1]] 

    v1_normal = dot([state1[3],state1[4]], normal)
    v2_normal = dot([state2[3],state2[4]], normal)
    v1_tangent = dot([state1[3],state1[4]], tangent)
    v2_tangent = dot([state2[3],state2[4]], tangent)
    
    v1_normal_after = ((ball1.mass - e[] *ball2.mass) * v1_normal + (1 + e[]) * ball2.mass *v2_normal) / (ball1.mass+ ball2.mass)
    v2_normal_after = ((ball2.mass - e[] * ball1.mass) * v2_normal+ (1 + e[])* ball1.mass * v1_normal) / (ball1.mass + ball2.mass)
    
    state1[3] = v1_normal_after *normal[1] + v1_tangent * tangent[1]
    state1[4] = v1_normal_after * normal[2] + v1_tangent* tangent[2]
    state2[3] = v2_normal_after *normal[1] + v2_tangent * tangent[1]
    state2[4] = v2_normal_after * normal[2] + v2_tangent* tangent[2]
end
#funkcja licząca kolizje ze zatrzymaną kulką
function collision_with_no_bouncing_ball(ball1::Ball,ball2::Ball,state::Vector{Float64})
    normal = normalize(ball2.position - ball1.position)
    velocity = [state[3],state[4]]
    v1_normal = dot(velocity, normal)
    velocity_new = velocity .- (1 + e[]) * v1_normal * normal 
    state[3],state[4] = velocity_new
end
#funkcja licząca kolizje ze ścianą
function collision_with_wall(ball1::Ball,state1::Vector{Float64}) 
    normal=normalize(ball1.position)
    v_normal = state1[3] * normal[1] + state1[4] * normal[2]
    if v_normal > 0 state1[3],state1[4] = [(state1[3]-(2) * v_normal * normal[1])*e[],(state1[4]-(2) * v_normal * normal[2])*e[]] end
end
#funkcja licząca deformacje ze ścianą
function deformation_with_wall(balls::Ball,time::Int,center::Observable,new_circle::Observable,x_circle::Observable,y_circle::Observable,initial_velocity::Float64)  
    alfa=min(norm(balls.velocity)/initial_velocity,2.0)
    deformation_per=k[]*alfa
    if deformation_per<0.02 return deformation_per end

    dist=norm(balls.position)
    if dist > 0.9
        normal=normalize(balls.position)
        if time==6
            balls.position -= (dist - 0.9)*normal
        elseif time>3
            balls.position -= 0.2*deformation_per*normal/3
        else
            balls.position += 0.2*deformation_per*normal/3
        end 
        center[]=[balls.position[1],balls.position[2]]
    end
    points = new_circle[]
    combined_points = [[points[1][i], points[2][i]] for i in 1:30] 
    temp_x=Vector{Float64}()
    temp_y=Vector{Float64}()
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
#funkcja licząca deformacje między kulami
function deformation_with_ball(ball1::Ball,ball2::Ball,time::Int,initial_velocity1::Float64,initial_velocity2::Float64,center1::Observable,center2::Observable)
    alfa1=min(norm(ball1.velocity-ball2.velocity)/initial_velocity1,2.0)
    alfa2=min(norm(ball1.velocity-ball2.velocity)/initial_velocity2,2.0)
    deformation_per1=(ball2.mass/(ball2.mass+ball1.mass))*k[]*alfa1
    deformation_per2=(ball1.mass/(ball2.mass+ball1.mass))*k[]*alfa2
    if deformation_per1+deformation_per2 <0.02 return deformation_per1+deformation_per2   end

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
#funkcja licząca korekte odległości kulek, aby nadmiarowo nie liczyło kolizji
function correction(state1::Vector{Float64},state2::Vector{Float64})
    dist=norm([state2[1],state2[2]] - [state1[1],state1[2]])
    normal=normalize([state2[1],state2[2]] - [state1[1],state1[2]])
    state1[1] -= (0.201-dist) * normal[1]/2
    state1[2] -= (0.201-dist) * normal[2]/2
    state2[1] += (0.201-dist) * normal[1]/2
    state2[2] += (0.201-dist) * normal[2]/2
end
#funkcja pełniąca role pętli symulacji
function start() 
    while true
        if !isopen(window.scene) || !is_simulating[] return end       
        for j in 1:length(balls)            
            #podawania wejścia do równania różniczkowego
            if no_bounces[j]
                input_vec[j].=[0.0,0.0]
                state_vec[j][3],state_vec[j][4]=[0.0,0.0] 
            else
                input_vec[j].=[0.0,-G*balls[j].mass]            
            end  
            #logika kolizji ze ścianą
            if check_collison(j) && !no_bounces[j]           
                no_deformation=false
                coll_time[j]+=1  

                if deformation_with_wall(balls[j],coll_time[j],center_vec[j],circle_points[j],x_circle[j],y_circle[j],initial_velocity[j]) < 0.02 no_deformation=true end 
                
                if coll_time[j]<6  && !no_deformation
                    continue
                elseif abs(state_vec[j][3]) <0.1 && abs(state_vec[j][4]) <0.2  && !no_bounces[j] && state_vec[j][2]<0
                    no_bounces[j]=true
                    balls[j].mass=9999999.0
                else
                    coll_time[j]=0
                    collision_with_wall(balls[j],state_vec[j])
                end  
            end
            #logika kolizji z kulkami
            for coll in j+1:length(balls)
                collision_occ=check_collision_balls(j,coll)
                
                if collision_occ 
                    no_deformation_ball=false
                    coll_time_ball[j]+=1 
                    coll_time_ball[coll]+=1 
                    if no_bounces[j] && no_bounces[coll] correction(state_vec[j],state_vec[coll]) end

                    if deformation_with_ball(balls[j],balls[coll],coll_time_ball[j],initial_velocity[j],initial_velocity[coll],center_vec[j],center_vec[coll]) < 0.02 no_deformation_ball=true end
                    
                    if  coll_time_ball[j]<6 && coll_time_ball[j]<6 && !no_deformation_ball continue end

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
            #liczenie równania różniczkowego
            if !check_collison(j) coll_time[j]=0 end
            if coll_time_ball[j]>6 coll_time_ball[j]=0 end
            if coll_time_ball[j]!=0
                continue
            else
                state_vec[j].=DiffEq(state_vec[j],input_vec[j],balls[j],STEP)  
            end
            #zabezpieczenie przed wylotem za okrąg
            if norm([state_vec[j][1],state_vec[j][2]])> 0.91
                dist=norm([state_vec[j][1],state_vec[j][2]])
                normal=normalize([state_vec[j][1],state_vec[j][2]])
                state_vec[j][1] -= (dist - 0.9101) *  normal[1]
                state_vec[j][2] -= (dist - 0.9101) *  normal[2]
            end
            #aktualizacja danych
            center_vec[j][] = [state_vec[j][1],state_vec[j][2]]
            balls[j].position=[state_vec[j][1],state_vec[j][2]] 
            balls[j].velocity=[state_vec[j][3],state_vec[j][4]]
        end 
        sleep(STEP)
    end             
end
#interfejs
window = Figure(size=(1920, 1080))
axis = Axis(window[1, 1],aspect=DataAspect())
spoint = select_point(axis.scene)
Makie.deactivate_interaction!(axis,:rectanglezoom)
lines!(axis,Circle(Point2f0(0.0,0.0), 1.0), color = :black, linewidth = 10)
menu_grid = GridLayout(window[1,2],tellwidth=false,tellheight=false)
mass_box = Textbox(menu_grid[1,1],placeholder="Mass",reset_on_defocus=true, height=30,width=150)
force_box = Textbox(menu_grid[2,1],placeholder="Force",reset_on_defocus=true,height=30,width=150)
submit = Button(menu_grid[3,1], label="Submit",height=30,width=150)
label_e = Label(menu_grid[4, 1], "e= $(e[])")
slider_e = Slider(menu_grid[5,1], range=0.01:0.01:1,startvalue=0.01)
label_k = Label(menu_grid[6, 1], "k= $(k[])")
slider_k = Slider(menu_grid[7,1], range=0:0.01:1,startvalue=0.0)
simulation = Button(menu_grid[8,1],label="Run Simulation",height=30,width=150)
reset_button = Button(menu_grid[9,1],label="Reset",height=30,width=150)
display(window)
#obsługa listenerów
#obsługa asynca oraz symulacji z pomocą Chata GPT oraz https://www.youtube.com/watch?v=L-gyDvhjzGQ&t=1s&ab_channel=JuliaDynamics
on(simulation.clicks) do simulate
    if length(balls)!=length(center_vec) return end
    if !is_simulating[]
        is_simulating[] = true
        simulation.label = "Simulating"
        @async begin
            if !isopen(window.scene) return end
            start()
            simulation.label = "Run Simulation" 
        end
    end  
end
on(reset_button.clicks) do reset
    is_simulating[]=false
    empty!(balls)
    empty!(state_vec)
    empty!(input_vec)
    empty!(no_bounces)
    empty!(initial_velocity)
    empty!(coll_time)
    empty!(coll_time_ball)
    empty!(center_vec)
    empty!(circle_points)
    empty!(x_circle)
    empty!(y_circle)
    create_state.control=1
    create_state.position_of_center=[0.0,0.0]
    create_state.alfa=0.0
    empty!(axis)
    lines!(axis,Circle(Point2f0(0.0,0.0), 1.0), color = :black, linewidth = 10)
end
on(spoint) do position 
    if !is_simulating[] click_plot_check(create_state,position) end
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
    if !is_simulating[]
        button_click_check(create_state, balls, mass_box, force_box)
        if length(balls)!=length(center_vec) && length(balls)>0
            create_state.control-=1
            pop!(balls)
            pop!(state_vec)
            pop!(input_vec)
            pop!(no_bounces)
            pop!(initial_velocity)
            pop!(coll_time)
            pop!(coll_time_ball)
        end
    end
end
