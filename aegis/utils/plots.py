#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Saturday March 23 14:29:53 2023

@authors: David Navarro, Antonio Santiago
"""

import colorlover as cl
import plotly.graph_objects as go

def hex_to_rgb(hex_string):
    hex_string = hex_string.lstrip("#")

    return tuple(int(hex_string[i:i + 2], 16) for i in (0, 2, 4))

palettes = {}
palettes["extreme"] = ["#8B0000", "#ADD8E6"]
palettes["purple"] = ["#F7B7A3", "#EA5F89", "#9B3192", "#57167E", "#2B0B3F"]
palettes["pastel"] = ["#A8D6B5", "#FFC9A5", "#C9CCE5", "#EFCEDF", "#DFF2B9"]

common_font = dict(family="Arial", size=18)

def pie_chart(labels, values, export_folder, tag, title, hovertext_labels:list=None, colours:str="purple"):
    colours = palettes[colours]
    fig = go.Figure(data=[go.Pie(labels=labels, sort=False, values=values, hoverinfo="label+value", textfont=common_font, hole=0.4,textinfo="percent")])
    fig.update_layout(title=title, title_font=dict(family="Arial", size=22), hoverlabel=dict(font_size=common_font["size"], font_family=common_font["family"]), title_x=0.5, legend=dict(font=common_font))
    if len(colours) == len(labels):
        pass
    elif len(colours) > len(labels):
        colours = colours[:len(labels)]
    else:
        rgb_colours = []
        for colour in colours:
            rgb_colours.append(hex_to_rgb(colour))
        colours = cl.interp(rgb_colours, len(labels))
    fig.update_traces(marker=dict(colors=colours, line=dict(color='#000000', width=0.05)))

    if hovertext_labels != None:
        fig.update_traces(hoverinfo='value+text', text=hovertext_labels)
        fig.update_layout(hoverlabel=dict(font_size=12))

    total = sum(values)
    fig.add_annotation(text=f'Total: {total}', x=0.5, y=0.5, font=common_font, showarrow=False)

    fig.write_image(f"{export_folder}{tag}.pdf")
    fig.write_html(f"{export_folder}{tag}.html")

def barplot(values, export_folder, tag, title, max_x:int=None):

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

    fig = go.Figure(data=[go.Histogram(x=values, marker_color='#4CAF50', marker_line_color='#000000', marker_line_width=0.05)])
    fig.update_layout(title=title, title_x=0.5, xaxis_title=tag, yaxis_title="Frequency", xaxis=dict(range=x_range), bargap=0.05, bargroupgap=0.1)
    fig.update_layout(title_font=dict(family="Arial", size=22), hoverlabel=dict(font_size=common_font["size"], font_family=common_font["family"]))
    fig.write_image(f"{export_folder}{tag}.pdf")
    fig.write_html(f"{export_folder}{tag}.html")
