import streamlit as st
import matplotlib
import matplotlib.pyplot as plt

from main1 import model_channels, plot_stuff

state = st.radio(
    "Sluice gate",
    ('Open', 'Closed'))

if (state =='Open'):
    s2 = 0
else:
    s2 = 1


sf, P_new, V = model_channels(s2)
fig = plot_stuff(sf, P_new, V)
st.pyplot(fig)