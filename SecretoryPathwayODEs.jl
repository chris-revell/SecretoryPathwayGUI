




# Define species concentrations
C_ER                # μM        ER-located collagen
C_H                 # μM        Collagen-HSP47 complex in ERGIC
C_G                 # μM        Golgi-located collagen
C_PG                # μM        Post-Golgi compartments-located collagen
C_PM                # μM        Plasma membrane-located collagen
C_E                 # μM        Extracellular collagen
H_ER                # μM        ER-located HSP47
H_G                 # μM        Golgi-located HSP47
K                   # μM        Pka
S                   # μM        SEC61A2
T                   # μM        TANGO1
P                   # μM        PDE4D
V                   # μM        VPS33B
M                   # μM        MMP14
D                   # μM        CTSK

# Define mass action kinetics parameters
μ_typ   = 0.0342    # μM        Median typical concentration
k₁      = 0.00822   # μM        SEC61A2 mean concentration
k₂      = 0.0105    # μM        TANGO1 mean concentration
k₃      = 0.0342    # μM        PDE4D mean concentration
k₆      = 0.00335   # μM        VPS33B mean concentration
k₇      = 0.00457   # μM        MMP14 mean concentration
k₁₅     = 0.0342    # μM        CTSK mean concentration
t̅       = 24.0      # h         Circadian period
k₄_a    = 4.84      # h⁻¹       HSP47 protein synthesis rate
k₄_b    = 0.00413   # h⁻¹       ER HSP47 degradation rate
k₄_c    = 0.00413   # h⁻¹       Golgi HSP47 degradation rate
k₅_a    = 1.89      # h⁻¹       PKA protein synthesis rate
k₅_b    = 0.0069    # h⁻¹       PKA protein degradation rate
k₈      = 0.0693    # h⁻¹       Col-I synthesis rate
k₁₀     = 3600.0    # h⁻¹       ERGIC to Golgi Col-I transition rate
k₁₂     = 0.0417    # h⁻¹       Post-golgi compartments to plasma membrane Col-I transition rate
k₅_c    = 29.2      # μM⁻¹      Rate constant for repression of PKA by PDE4D
k₋₉     = 1.22      # μM⁻¹h⁻¹   Col-I-HSP47 dissociation rate
k₁₁     = 1.22      # μM⁻¹h⁻¹   Golgi to post-golgi compartments Col-I transition rate
k₁₃     = 12.2      # μM⁻¹h⁻¹   Plasma membrane to extracellular Col-I transition rate
k₁₄     = 1.22      # μM⁻¹h⁻¹   CTSK-dependent extracellular Col-I degradation rate
k₁₆     = 1.22      # μM⁻¹h⁻¹   Golgi to ER HSP47 transition rate
k₉      = 0.356     # μM⁻²h⁻¹   Col-I-HSP47 complex association rate
G̅                   # Unitless  Mean pulse function value determined by mean value theorem
t₀                  # h         Phase
