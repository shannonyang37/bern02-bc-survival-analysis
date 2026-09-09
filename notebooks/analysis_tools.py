# %% [markdown]
# ### Notebook 2: Defining Analysis Tools
# 
# This notebook serves as the backend graphical and statistical engine for the clinical analysis. It contains reusable, modular functions designed to automate univariate exploratory data analysis (EDA) for time-to-event clinical data.
# 
# It contains:
# * Automated statistical routing (Welch's T-Test, ANOVA, Chi-Square, Log-Rank).
# * Generation of 3-panel clinical dashboards (Absolute Outcomes, Mortality Velocity, Kaplan-Meier Trajectories).
# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import stats
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
# %%
def plot_box_numerical_with_stats(column_name, dataframe):
    """
    Generates a boxplot for a given numerical clinical factor, coupled with a Welch's T-Test (two-sample independent t-test) to validate the result.

    Parameters:
    ------------
    column_name: str
        The exact string name of the column in the dataframe to be analyzed (e.g. 'Age', 'Tumor Size')
    dataframe: pandas.DataFrame
        The cleaned clinical dataset containing the numerical columns and a 'Status' column

    Outputs:
    --------
    Displays a matplotlib figure of the boxplot with statistical annotations (p-value).
    """

    # Isolate the numerical arrays for the two groups ('Alive' & 'Dead')
    alive = dataframe[dataframe['Status'] == 'Alive'][column_name]
    dead = dataframe[dataframe['Status'] == 'Dead'][column_name]

    # Carry out the Welch's T-Test to calculate the p-value
    t_stat, p_val = stats.ttest_ind(alive, dead, equal_var=False)

    # Format the p-value to look clean (e.g., "p < 0.0001" or "p = 0.0312")
    if p_val < 0.0001:
        p_text = "p < 0.0001"
    else:
        p_text = f"p = {p_val:.4f}"

    # Draw the boxplot
    plt.figure(figsize=(5, 6))
    sns.boxplot(x='Status', y=column_name, data=dataframe, palette='pastel', hue='Status')

    # Styling and Layout
    plt.suptitle(f'Distribution of {column_name} by Survival', fontsize=12, fontweight='bold', y=0.96)
    plt.title(f'(Welch\'s T-Test: {p_text})', fontsize=10, pad=10)
    plt.xlabel('Patient Survival Status', fontsize=10)
    plt.ylabel(column_name, fontsize=10)

    plt.tight_layout()
    plt.show()
# %%
def survival_ttest(column_name, order, dataframe):
    """
    Carries out the Welch's T-Test to validate the mortality velocity for deceased patients across a categorical clinical factors containing exactly two groups.

    Parameters:
    -----------
    column_name : str
        The exact string name of the categorical column to be analyzed (e.g., 'Estrogen Status').
    order : list of str
        The exact two categories to be tested, in the desired order (e.g., ['Positive', 'Negative']).
    dataframe : pandas.DataFrame
        The cleaned clinical dataset containing the categorical column and 'Survival Months'.

    Outputs:
    --------
    t_stat: float
        The calculated Welch's T-statistic.
    p_val: float
        The calculated p-value.
    """

    # Safety check: Welch's Test requires exactly two groups to run
    if len(order) != 2:
        raise ValueError(f"Welch's Test requires exactly two groups. '{column_name}' has {len(order)}.")

    # Isolate the numerical arrays for the two groups
    group_one = dataframe[dataframe[column_name] == order[0]]['Survival Months']
    group_two = dataframe[dataframe[column_name] == order[1]]['Survival Months']

    # Carry out the Welch's T-Test to calculate the p-value
    t_stat, p_val = stats.ttest_ind(group_one, group_two, equal_var=False)

    return t_stat, p_val
# %%
def survival_anova(column_name, order, dataframe, numerical_column='Survival Months'):
    """
    Carries out a One-Way ANOVA to validate the mortality velocity for deceased patients across a categorical clinical factor containing 3 or more groups.

    Parameters:
    -----------
    column_name : str
        The exact string name of the categorical column to be analyzed (e.g., 'Grade').
    order : list of str
        The exact categories to be tested, in the desired order (e.g., ['Grade I', 'Grade II', 'Grade III']).
    dataframe : pandas.DataFrame
        The cleaned clinical dataset.
    numerical_column : str, default='Survival Months'
        The continuous variable to be tested across the groups.

    Returns:
    --------
    f_stat : float
        The calculated F-statistic.
    p_val : float
        The calculated p-value.
    """

    # Extract groups based explicitly on the provided order
    raw_groups = [dataframe[dataframe[column_name] == category][numerical_column].dropna() for category in order]

    # Safety Check: Filter out any completely empty groups to prevent crashes
    valid_groups = [group for group in raw_groups if not group.empty]

    # Safety Check: ANOVA requires at least 2 populated groups to compare
    if len(valid_groups) < 2:
        raise ValueError(f"ANOVA requires at least 2 populated groups. '{column_name}' only has {len(valid_groups)}.")

    # 4. Run the ANOVA using the '*' unpacking operator
    f_stat, p_val = stats.f_oneway(*valid_groups)

    return f_stat, p_val
# %%
def survival_logrank(column_name, order, dataframe):
    """
    Carries out a Multivariate Log-Rank test to validate the statistical significance of Kaplan-Meier survival trajectories across specified categorical groups.

    Parameters:
    -----------
    column_name : str
        The exact string name of the categorical column to be analyzed (e.g., 'T Stage').
    order : list of str
        The exact categories to be tested (e.g., ['T1', 'T2', 'T3', 'T4']).
    dataframe : pandas.DataFrame
        The cleaned clinical dataset containing 'Status' and 'Survival Months'.

    Returns:
    --------
    p_val : float
        The calculated Log-Rank p-value.
    """

    # Create a lightweight temporary copy containing only what we need
    cols_needed = [column_name, 'Status', 'Survival Months']
    df_temp = dataframe[cols_needed].copy()

    # Create the binary Event column (1=Dead, 0=Alive)
    df_temp['_Event'] = df_temp['Status'].map({'Alive': 0, 'Dead': 1})

    # Filter the data based on the order list
    df_temp = df_temp[df_temp[column_name].isin(order)]

    # Drop missing values to prevent crashes
    df_temp = df_temp.dropna(subset=[column_name, 'Survival Months', '_Event'])

    # Safety Check: Ensure we have at least 2 populated groups left to compare
    if df_temp[column_name].nunique() < 2:
        raise ValueError(f"Log-Rank test requires at least 2 populated groups. '{column_name}' lacks sufficient data.")

    # Run the Multivariate Log-Rank Test
    results = multivariate_logrank_test(
        event_durations=df_temp['Survival Months'],
        groups=df_temp[column_name],
        event_observed=df_temp['_Event']
    )

    return results.p_value
# %%
def plot_comprehensive_analysis_dashboard(column_name, order, dataframe):
    """
    Generates a 3-panel clinical dashboard, coupled with statistical tests to validate the results.

    Parameters:
    -----------
    column_name : str
        The exact string name of the categorical column to be analyzed (e.g., 'T Stage').
    order : list of str
        The exact categories to be tested (e.g., ['T1', 'T2', 'T3', 'T4']).
    dataframe : pandas.DataFrame
        The cleaned clinical dataset containing 'Status' and 'Survival Months'.

    Outputs:
    --------
    Top Left: Bar chart (Absolute Survival & Mortality Rates with Chi-Square)
    Top Right: Boxplot (Timeline to Mortality for Deceased with ANOVA/T-Test)
    Bottom: Kaplan-Meier Survival Curve (Long-term probability with Log-Rank Test)
    """

    # Setup the Centered Layout Canvas
    fig = plt.figure(figsize=(16, 11))

    # Create a 2-row, 4-column grid
    gs = gridspec.GridSpec(2, 4, height_ratios=[1, 1])

    # Ax1: Top Left (Row 0, spans columns 0 to 2)
    ax1 = fig.add_subplot(gs[0, 0:2])

    # Ax2: Top Right (Row 0, spans columns 2 to 4)
    ax2 = fig.add_subplot(gs[0, 2:4])

    # Ax3: Bottom Center (Row 1, spans columns 1 to 3)
    ax3 = fig.add_subplot(gs[1, 1:3])

    # Create a local temporary event column
    df_temp = dataframe.copy()
    if '_Event' not in df_temp.columns:
        df_temp['_Event'] = df_temp['Status'].map({'Alive': 0, 'Dead': 1})

    # ==========================================
    # PLOT 1 (Top Left): Absolute Patient Outcomes
    # ==========================================
    # Chi-Square Test
    contingency_table = pd.crosstab(df_temp[column_name], df_temp['Status'])
    chi2_stat, p_val_chi, dof, expected_frequencies = stats.chi2_contingency(contingency_table)
    p_text_chi = "p < 0.0001" if p_val_chi < 0.0001 else f"p = {p_val_chi:.4f}"

    # Draw Bar Chart
    sns.countplot(x=column_name, data=df_temp, palette='pastel', hue='Status', order=order, ax=ax1)

    # Add Percentage Labels
    category_totals = {cat: df_temp[df_temp[column_name] == cat].shape[0] for cat in order}
    for container in ax1.containers:
        labels = []
        for idx, bar in enumerate(container):
            height = bar.get_height()
            current_category = order[idx]
            total = category_totals.get(current_category, 0)
            rate = (height / total) * 100 if total > 0 else 0
            labels.append(f'{int(height)}\n({rate:.1f}%)' if height > 0 else '')
        ax1.bar_label(container, labels=labels, padding=5, fontsize=8.5, fontweight='bold')

    # Styling and layout
    ax1.set_title(f'Survival & Mortality Rates\n[Chi-Square: {p_text_chi}]', fontsize=12, fontweight='bold', pad=10)
    ax1.set_xlabel(column_name, fontsize=10)
    ax1.set_ylabel('Patient Count', fontsize=10)
    ax1.set_ylim(0, df_temp[column_name].value_counts().max() * 1.15)
    ax1.grid(True, axis='y', alpha=0.3, linestyle='--')

    # Rotation for long groups names
    if max([len(str(label)) for label in order]) > 10:
        ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha='right')

    # ==========================================
    # PLOT 2 (Top Right): Velocity of Mortality
    # ==========================================
    # Generate a dataframe with only deceased patients
    deceased_df = df_temp[df_temp['Status'] == 'Dead']

    # Statistical Test: Welch's Test or ANOVA according to number of groups
    if len(order) == 2:
        t_stat, p_val_box = survival_ttest(column_name, order, deceased_df)
        test_name_box = "Welch's T-Test"
    else:
        f_stat, p_val_box = survival_anova(column_name, order, deceased_df)
        test_name_box = "ANOVA"

    p_text_box = f"{test_name_box}: p < 0.0001" if p_val_box < 0.0001 else f"{test_name_box}: p = {p_val_box:.4f}"

    # Draw Boxplot
    sns.boxplot(x=column_name, y='Survival Months', data=deceased_df, order=order, palette='pastel', ax=ax2)

    # Styling and layout
    ax2.set_title(f'Mortality Velocity (Deceased Cohort)\n[{p_text_box}]', fontsize=12, fontweight='bold', pad=10)
    ax2.set_xlabel(column_name, fontsize=10)
    ax2.set_ylabel('Months Survived', fontsize=10)
    ax2.grid(True, axis='y', alpha=0.3, linestyle='--')

    # Rotation for long groups names
    if max([len(str(label)) for label in order]) > 10:
        ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45, ha='right')

    # ==========================================
    # PLOT 3 (Bottom): Kaplan-Meier Curve
    # ==========================================
    # Set a color for the plot
    custom_colors = sns.color_palette('pastel', n_colors=len(order))

    # Draw Lines
    for i, category in enumerate(order):
        category_group = df_temp[df_temp[column_name] == category]
        if category_group.empty:
            continue

        kmf = KaplanMeierFitter()
        kmf.fit(category_group['Survival Months'], category_group['_Event'], label=str(category))
        kmf.plot_survival_function(ax=ax3, color=custom_colors[i], ci_show=False, linewidth=2.5)

    # Statistical Test: Log-Rank Test
    logrank_p_val = survival_logrank(column_name, order, dataframe)
    p_text_km = "Log Rank Test: p < 0.0001" if logrank_p_val < 0.0001 else f"Log Rank Test: p = {logrank_p_val:.4f}"

    # Styling and layout
    ax3.set_title(f'Long-Term Survival Trajectories (Kaplan-Meier)\n[{p_text_km}]', fontsize=12, fontweight='bold', pad=10)
    ax3.set_xlabel('Months Since Diagnosis', fontsize=10)
    ax3.set_ylabel('Probability of Survival (0.0 to 1.0)', fontsize=10)
    ax3.set_ylim(0, 1.05)
    ax3.grid(True, alpha=0.3, linestyle='--')

    # ==========================================
    # Final Layout Adjustments
    # ==========================================
    plt.suptitle(f'Comprehensive Clinical Analysis: {column_name}', fontsize=20, fontweight='bold', y=1.02)

    plt.tight_layout()
    plt.show()
# %% [markdown]
# Save this notebook as a Python file to use the functions in future notebooks.
