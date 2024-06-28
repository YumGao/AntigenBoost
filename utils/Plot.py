import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

def Dot_plot(df, filename='plot.png'):
    # Reset index to start from 0
    df.reset_index(drop=True, inplace=True)
    
    # Create the figure and assign it to a variable
    fig = plt.figure(figsize=(8, 9))
    
    # Scatter plot for all strategies
    plt.scatter(df['Conservative_Score'], 
                df['ΔSeqScore'], 
                s=50, c='deepskyblue', label='Strategies')
    
    # Adding labels     
    for i, txt in enumerate(df["strategy"]):
        plt.annotate(txt, (df["Conservative_Score"][i], df["ΔSeqScore"][i]), fontsize=10, ha='right')
    
    # Labels and title
    plt.xlabel("Conservative Score")
    plt.ylabel("ΔSeqScore")
    #plt.title("Dot Plot")
    
    # Adjusting the limits of the axes to leave more space at the edge
    #plt.xlim(df['Conservative_Score'].min() - 1, df['Conservative_Score'].max() + 1)
    
    # Set transparency and background color
    fig.patch.set_alpha(0.0)
    plt.gca().set_facecolor('none')
    
    # Set the y-axis grid with a step of 0.2
    plt.gca().yaxis.set_major_locator(MultipleLocator(0.2))

    plt.legend(loc='upper right', bbox_to_anchor=(1.52, 1))

    plt.grid(True)

    # Save the figure to a PNG file
    plt.savefig(filename, format='png', bbox_inches='tight')
    print(f"Saved score analysis plot into file {filename}")

# Example usage (assuming you have a DataFrame `df`):
# df = pd.DataFrame({
#     'strategy': ['K320Q', 'K320F', 'K320M', 'Other1', 'Other2'],
#     'Conservative_Score': [1.0, 2.0, 3.0, 4.0, 5.0],
#     'ΔSeqScore': [0.5, -0.3, 0.1, 0.4, 0.2]
# })
# Dot_plot(df, 'my_plot.png')


