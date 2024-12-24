import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from scipy.stats import linregress, pearsonr, spearmanr

# Load the updated data
updated_data = pd.read_csv("linear.csv", encoding="latin1")  # Replace with your file path

# Clean the updated data
clean_updated_data = updated_data[['hydrophobicity (H)', 'b-sheet propensity (B)', 'charge (C)', 'UL (U)']].dropna()

# Define predictors (X) and response variable (y)
X_updated = clean_updated_data[['hydrophobicity (H)', 'b-sheet propensity (B)', 'charge (C)']].values
y_updated = clean_updated_data['UL (U)'].values

# Add an intercept column for linear regression
X_updated_intercept = np.column_stack((X_updated, np.ones(X_updated.shape[0])))

# Perform transformations: original, log, square root
y_updated_log = np.log10(y_updated)
y_updated_sqrt = np.sqrt(y_updated)


# Function to calculate R^2, regression coefficients, slope, r, and p-values
def analyze_model(X, y):
    # Train the regression model
    model = LinearRegression(fit_intercept=False)
    model.fit(X, y)
    y_pred = model.predict(X)
    
    # Get regression coefficients
    coefficients = model.coef_
    
    # Calculate R^2
    r2 = r2_score(y, y_pred)
    
    # Calculate correlation coefficient and p-value for observed vs predicted
    r_value, p_value = pearsonr(y, y_pred)
    
    # Calculate the slope of the line of best fit in predicted vs observed space
    slope, intercept, _, _, _ = linregress(y, y_pred)
    
    return r2, slope, r_value, p_value, y_pred, coefficients

# Analyze each transformation
results = {}
for label, y in [("Linear", y_updated), ("Log", y_updated_log), ("Square Root", y_updated_sqrt)]:
    r2, slope, r, p, y_pred, coefficients = analyze_model(X_updated_intercept, y)
    results[label] = {
        "R^2": r2,
        "Slope": slope,
        "r": r,
        "p-value": p,
        "Predicted": y_pred,
        "Coefficients": coefficients
    }

# Generate publication-grade plots for each model
def plot_with_annotations(observed, predicted, model_name, slope, r_value, p_value, r2_value):
    plt.figure(figsize=(8, 7))  # Square plot
    plt.scatter(observed, predicted, color='black', s=30, label="Data Points")  # Dots for data points
    
    # Calculate and plot the regression line
    x_vals = np.linspace(observed.min(), observed.max(), 100)
    regression_line = slope * x_vals + np.mean(predicted - slope * observed)
    plt.plot(x_vals, regression_line, color="blue", linestyle="-", linewidth=1.5, label="Regression Line")
    
    # Add annotations
    plt.text(
        0.05, 0.8, 
        f"$R^2 = {r2_value:.3f}$\n$r = {r_value:.3f}$\n$p = {p_value:.4g}$\nSlope = {slope:.3f}", 
        transform=plt.gca().transAxes, fontsize=15, verticalalignment='top',
        bbox=dict(boxstyle="square", edgecolor="black", facecolor="white")
    )
    
    # Calculate shared axis limits with margins
    min_val = min(observed.min(), predicted.min())
    max_val = max(observed.max(), predicted.max())
    margin = (max_val - min_val) * 0.05  # 5% margin on each side
    plt.xlim(min_val - margin, max_val + margin)
    plt.ylim(min_val - margin, max_val + margin)
    
    # Ensure equal aspect ratio
    plt.gca().set_aspect('equal', adjustable='box')  # Ensures axes are scaled equally
    
    # Add the perfect fit line
    #plt.plot([min_val, max_val], [min_val, max_val], color="red", linestyle="--", linewidth=1.5, label="Perfect Fit Line")
    
    # Customize axis labels and ticks
    plt.xlabel(f"Observed Values ({model_name})", fontsize=17)
    plt.ylabel(f"Predicted Values ({model_name})", fontsize=17)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    
    # Customize grid and legend
    plt.grid(False)  # No grid for cleaner look
    plt.legend(frameon=False, fontsize=12, loc="upper left")  # Updated legend with descriptions
    plt.tight_layout()  # Optimize spacing
    
    # Save and display the plot
    plt.savefig('regression_with_margin.svg', dpi=600, bbox_inches='tight')
    plt.show()


# Generate plots and print results
for label, result in results.items():
    observed = y_updated if label == "Linear" else (y_updated_log if label == "Log" else y_updated_sqrt)
    predicted = result["Predicted"]
    slope = result["Slope"]
    r_value = round(result["r"],2)
    p_value = round(result["p-value"],4)
    r2_value = round(result["R^2"],2)
    coefficients = result["Coefficients"]
    
    print(f"\n{label} Model Coefficients:")
    print(f"  Coefficients: {coefficients}")
    
    plot_with_annotations(observed, predicted, label, slope, r_value, p_value, r2_value)

# Print summary of results
print("Regression Analysis Results:")
for model_name, metrics in results.items():
    print(f"\n{model_name} Model:")
    print(f"  R^2: {metrics['R^2']:.4f}")
    #print(f"  Slope: {metrics['Slope']:.4f}")
    print(f"  Correlation Coefficient (r): {metrics['r']:.4f}")
    print(f"  P-Value: {metrics['p-value']:.4g}")
    print(f"  Coefficients: {metrics['Coefficients']}")
