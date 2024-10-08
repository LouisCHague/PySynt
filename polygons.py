import matplotlib.pyplot as plt
import matplotlib.patches as patches

def draw_segments_with_fill(RStart, Rend, Qstart, Qend, R_y=8, Q_y=4):
    # Create a figure and axis
    fig, ax = plt.subplots()
    
    # Define the corners of the filled area
    top_left = (RStart, R_y)
    top_right = (Rend, R_y)
    bottom_left = (Qstart, Q_y)
    bottom_right = (Qend, Q_y)
    
    # Create a polygon to fill the area between the lines
    filled_area = patches.Polygon([top_left, top_right, bottom_right, bottom_left], 
                                  closed=True, facecolor='lightblue', edgecolor='none', alpha=0.5)
    
    # Add the threads to the plot
    ax.add_patch(filled_area)
    
    ax.set_xlim(min(RStart, Qstart) - 1, max(Rend, Qend) + 1)
    ax.set_ylim(Q_y - 2, R_y + 2)
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    ax.set_title('Filled Area Between Segments Without Top and Bottom Lines')
    
    # Add grid and reference lines
    plt.grid()
    plt.axhline(0, color='black', linewidth=0.5, ls='--')
    plt.axvline(0, color='black', linewidth=0.5, ls='--')
    
    # Display the plot
    plt.show()

# Example usage
draw_segments_with_fill(RStart=1, Rend=2, Qstart=2, Qend=3)