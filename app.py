import numpy as np
import plotly.graph_objects as go

from viktor import ViktorController
from viktor.parametrization import ViktorParametrization, NumberField, Text
from viktor import ViktorController
from viktor.geometry import SquareBeam
from viktor.views import GeometryView, GeometryResult
from viktor.views import PlotlyView, PlotlyResult
from pathlib import Path
from viktor.views import PDFView, PDFResult

def I_3(a, b, z):

  alfa_1 = np.arctan((a+b)/z) - np.arctan(b/z)
  alfa_2 = np.arctan(b/z)
  result = (1/np.pi) * (((a+b)/a)*(alfa_1+alfa_2) - (b/a)*alfa_2)
  
  return result

def tramo_1(d, s, h, z, gamma, d_1):

  I_3_1 = I_3(s*h, d_1, z)
  I_3_2 = I_3(s*h, d-d_1, z)
  q = gamma * h
  s_total = q * (I_3_1 + I_3_2)

  return s_total


def tramo_2(d, s, h, z, gamma, d_2):
  
  I_3_1 = I_3(s*h-d_2, 0, z)
  q_1 = ((s*h-d_2)/s) * gamma
  I_3_2 = I_3(s*h, d+d_2, z)
  q_2 = h * gamma
  I_3_3 = I_3(d_2, 0, z)
  q_3 = (d_2/s) * gamma
  s_total = I_3_1*q_1 + I_3_2*q_2 - I_3_3*q_3

  return s_total

def tramo_3(d, s, h, z, gamma, d_3):
    
  I_3_1 = I_3(s*h, d + s*h + d_3, z)
  I_3_2 = I_3(s*h, d_3, z)
  q = gamma * h
  s_total = q * (I_3_1-I_3_2)
  
  return s_total

def stress_line(d, s, h, z, gamma):
    
  D_1 = np.arange(0, d/2 + 1e-5, 0.25)
  S_1 = tramo_1(d, s, h, z, gamma, D_1)

  D_2 = np.arange(0, h*s + 1e-5, 0.25)
  S_2 = tramo_2(d, s, h, z, gamma, D_2)

  D_3 = np.arange(0, h*2 + 1e-5, 0.25)
  S_3 = tramo_3(d, s, h, z, gamma, D_3)

  D = np.concatenate((-D_1[::-1] , D_2[1:-1], D_3 + h*s), axis=0) + d/2
  S = np.concatenate((S_1[::-1], S_2[1:-1], S_3), axis=0)

  DD = np.concatenate((-D[::-1] , D[1:]), axis=0)
  SS = np.concatenate((S[::-1] , S[1:]), axis=0)

  return DD, SS

def stress_mesh(d, s, h, z_max, gamma):
    
  y_values = np.arange(0.1, z_max + 1e-5, 0.25)
  S_values = []
  
  for element in y_values:
    S_values.append(stress_line(d, s, h, element, gamma)[1])
    x_values = stress_line(d, s, h, element, gamma)[0]
    
  return x_values, y_values, S_values

class Parametrization(ViktorParametrization):
    text_01 = Text("# Vertical stress in a semi-infinite mass due to enbankment loading")
    text_02 = Text("This program calculate the vertical stresses of a dam and a 2D colormap of this")
    
    width = NumberField("Crest Width", default=16, step=0.1, suffix="m")
    slope = NumberField("Slope", default=3, step=0.1, suffix="")
    height = NumberField("Height", default=50, step=0.1, suffix="m")
    depth = NumberField("Depth", default=70, step=0.1, suffix="m")
    gamma = NumberField("γ", default=19, step=0.1, suffix="kN/m3")
    
class Controller(ViktorController):
    viktor_enforce_field_constraints = True  # Resolves upgrade instruction https://docs.viktor.ai/sdk/upgrades#U83

    label = 'My Beam'
    parametrization = Parametrization
    
    @PlotlyView("My Graph", duration_guess=1)
    def stress_plot(self, params, **kwargs):
        # Results of the mesh
        results = stress_mesh(params.width, 
                              params.slope, 
                              params.height, 
                              params.depth, 
                              params.gamma)     
        # Plot
        fig = go.Figure()
        fig.add_trace(
            go.Contour(
                z=results[2],
                x=results[0],
                y=results[1]
            )
        )
        
        d = params.width
        s = params.slope
        h = params.height
        gamma = params.gamma
        
        fig.add_trace(
            go.Scatter(
                x=[-d/2 - h*s, d/2 + h*s, d/2, -d/2, -d/2 - h*s], 
                y=[0, 0, -h, -h, 0],
                fill="toself",
                fillcolor="rgba(201,201,201,1)",
                mode='lines',
                line=dict(color="gray", width=2)
                )
        )

        fig.add_annotation(
                x= d/2 + h/s + 1.5 * h,
                y= - h,
                xref="x",
                yref="y",
                text=
                "Height = "+ str(round(float(h), 1)) +
                "m <br>"+
                "Crest Width = " + str(round(float(d), 1)) +
                "m <br>"+ 
                "Slope = "+ str(round(float(s), 1)) +
                "<br>"+ 
                "U.W. = " + str(round(float(gamma), 1)) +
                "kN/m3",
                showarrow=False,
                font=dict(
                    family="Courier New, monospace",
                    size=11,
                    color="gray"
                    ),
                align="center",
                bordercolor="gray",
                borderwidth=2,
                borderpad=4,
                bgcolor="rgba(220,220,220,1)",
                opacity=0.8
                )

        fig.update_layout(yaxis= {'autorange': 'reversed'})
        fig.update_layout(xaxis_title = "Distance (m)",
                        yaxis_title = "Depth (m)",
                        title = "Vertical stress in a semi-infinite mass due to enbankment loading")

        return PlotlyResult(fig.to_json())   

    @PDFView("Theory", duration_guess=1)
    def get_pdf_view(self, params, **kwargs):
        print(Path(__file__).parent)
        file_path = Path(__file__).parent / 'doc.pdf'
        return PDFResult.from_path(file_path)       