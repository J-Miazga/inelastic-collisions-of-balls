using GeometryBasics
using GLMakie

# Function to generate the coordinates for a circle
function generate_circle_points(center, radius, n_points = 100)
    angles = LinRange(0, 2pi, n_points)  # Create angle steps
    xs = center[1] .+ radius .* cos.(angles)  # X coordinates
    ys = center[2] .+ radius .* sin.(angles)  # Y coordinates
    return xs, ys
end

ox=Observable(2)
oy=Observable(2)
# Create a figure and axis
fig = Figure(size = (1920, 1080))
ax = Axis(fig[1, 1], aspect= DataAspect())

# Define the circle parameters
center = Point2f0(0.0, 0.0)  # Circle center
radius = 1.0  # Circle radius

# Generate circle points
xs, ys = generate_circle_points(center, radius)

# Plot the circle using lines
#lines!(ax, xs, ys, color = :black, linewidth = ox)
lines!(ax, Circle(Point2f0(0.0, 0.0), 1.0), color = :black, linewidth = 10)
scatter!(ax,ox, oy,color=:black)

# Show the figure
display(fig)