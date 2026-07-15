#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Saturday March 23 14:29:53 2023

@authors: David Navarro, Antonio Santiago
"""

import colorlover as cl
import plotly.graph_objects as go
import matplotlib.pyplot as plt

def hex_to_rgb(hex_string):
    hex_string = hex_string.lstrip("#")
    return tuple(int(hex_string[i:i + 2], 16) for i in (0, 2, 4))

palettes = {}
palettes["extreme"] = ["#8B0000", "#ADD8E6"]
palettes["purple"] = ["#F7B7A3", "#EA5F89", "#9B3192", "#57167E", "#2B0B3F"]
palettes["pastel"] = ["#A8D6B5", "#FFC9A5", "#C9CCE5", "#EFCEDF", "#DFF2B9"]

common_font = dict(family="Arial", size=18)

def parse_to_matplotlib_colors(colours):
    """Converts hex or rgb() strings from colorlover into a format matplotlib understands."""
    matplotlib_colors = []
    for c in colours:
        if isinstance(c, str) and c.startswith("rgb"):
            # Parse 'rgb(r,g,b)' into normalized float tuples (r/255, g/255, b/255)
            parts = c.replace("rgb(", "").replace(")", "").split(",")
            matplotlib_colors.append(tuple(float(p.strip()) / 255.0 for p in parts))
        else:
            matplotlib_colors.append(c)
    return matplotlib_colors

def save_pdf_pie(labels, values, colours, export_folder, tag, title):
    """Generates the PDF pie chart using Matplotlib."""
    matplotlib_colors = parse_to_matplotlib_colors(colours)
    total = sum(values)

    fig, ax = plt.subplots(figsize=(6, 6))
    
    # Create the donut chart (wedge width < 1 creates the center hole)
    wedges, texts, autotexts = ax.pie( #type: ignore
        values,
        labels=labels,
        colors=matplotlib_colors,
        autopct='%1.1f%%',
        startangle=90,
        counterclock=False,
        wedgeprops=dict(width=0.6, edgecolor='#000000', linewidth=0.5),
        textprops=dict(family="sans-serif", size=12)
    )
    
    # Add the total value inside the center hole
    ax.text(0, 0, f'Total:\n{total}', ha='center', va='center', fontsize=14, family="sans-serif", weight='bold')
    
    ax.set_title(title, fontsize=16, family="sans-serif", pad=20)
    plt.tight_layout()
    plt.savefig(f"{export_folder}{tag}.pdf", format='pdf', bbox_inches='tight')
    plt.close(fig)

def save_pdf_barplot(values, min_x, max_x, export_folder, tag, title):
    """Generates the PDF histogram using Matplotlib."""
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Matplotlib's hist function
    ax.hist(values, bins='auto', range=(min_x, max_x), color='#4CAF50', edgecolor='#000000', linewidth=0.5, rwidth=0.95)
    
    ax.set_title(title, fontsize=16, family="sans-serif", pad=15)
    ax.set_xlabel(tag, fontsize=12, family="sans-serif")
    ax.set_ylabel("Frequency", fontsize=12, family="sans-serif")
    ax.set_xlim(min_x, max_x)
    
    # Simple grid for readability
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    ax.set_axisbelow(True)
    
    plt.tight_layout()
    plt.savefig(f"{export_folder}{tag}.pdf", format='pdf', bbox_inches='tight')
    plt.close(fig)


def pie_chart(labels:list[str], values:list[int], export_folder:str, tag:str, title:str, hovertext_labels:list|None=None, palette_name:str="purple"):
    colours = palettes[palette_name]
    
    if len(colours) == len(labels):
        pass
    elif len(colours) > len(labels):
        colours = colours[:len(labels)]
    else:
        rgb_colours = []
        for colour in colours:
            rgb_colours.append(hex_to_rgb(colour))
        colours = cl.interp(rgb_colours, len(labels))

    # 1. Generate Interactive HTML using Plotly
    fig = go.Figure(data=[go.Pie(labels=labels, sort=False, values=values, hoverinfo="label+value", textfont=common_font, hole=0.4, textinfo="percent")])
    fig.update_layout(title=title, title_font=dict(family="Arial", size=22), hoverlabel=dict(font_size=common_font["size"], font_family=common_font["family"]), title_x=0.5, legend=dict(font=common_font))
    fig.update_traces(marker=dict(colors=colours, line=dict(color='#000000', width=0.05)))

    if hovertext_labels is not None:
        fig.update_traces(hoverinfo='value+text', text=hovertext_labels)
        fig.update_layout(hoverlabel=dict(font_size=12))

    total = sum(values)
    fig.add_annotation(text=f'Total: {total}', x=0.5, y=0.5, font=common_font, showarrow=False)
    fig.write_html(f"{export_folder}{tag}.html")

    # 2. Generate Static PDF using Matplotlib (No Chrome/Kaleido dependency)
    save_pdf_pie(labels, values, colours, export_folder, tag, title)


def barplot(values:list[int], export_folder:str, tag:str, title:str, max_x:int|None=None):
    if max_x == None:
        max_x = max(values) + 1
    else:
        tag += f"_max_x_{max_x}"

    min_x = 0
    if len(values) > 0:
        if min(values) >= 0:
            min_x = 0
        else:
            min_x = min(values) - 1
    x_range = [min_x, max_x]

    # 1. Generate Interactive HTML using Plotly
    fig = go.Figure(data=[go.Histogram(x=values, marker_color='#4CAF50', marker_line_color='#000000', marker_line_width=0.05)])
    fig.update_layout(title=title, title_x=0.5, xaxis_title=tag, yaxis_title="Frequency", xaxis=dict(range=x_range), bargap=0.05, bargroupgap=0.1)
    fig.update_layout(title_font=dict(family="Arial", size=22), hoverlabel=dict(font_size=common_font["size"], font_family=common_font["family"]))
    fig.write_html(f"{export_folder}{tag}.html")

    # 2. Generate Static PDF using Matplotlib (No Chrome/Kaleido dependency)
    save_pdf_barplot(values, min_x, max_x, export_folder, tag, title)