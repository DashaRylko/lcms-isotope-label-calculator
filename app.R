
library(shiny)
library(tidyverse)
library(bslib) 

# Isotopes dictionary
izotopes_1 = list(
  "углерод"= list(`12`= 0.9894, `13.003354835`= 0.0106),
  "водород"= list(`1.0078250322`= 0.999855, `2.0141017781`= 0.000145),
  "азот"= list(`14.003074004`= 0.996205, `15.000108899`= 0.003795),
  "кислород"= list(`15.994914619`= 0.99757, `16.999131757`= 0.0003835, `17.999159613`= 0.002045),
  "сера"= list(`31.972071174`= 0.9485, `32.97145891`= 0.00763, `33.9678670`= 0.04365, `35.967081`= 0.000158),
  "хлор"= list(`34.9688527`= 0.758, `36.9659026`= 0.242),
  "бром"= list(`78.918338`= 0.5065, `80.916288`= 0.4935),
  "фтор"=list(`18.998403162`=1),
  "фосфор"=list(`30.973761998`=1)
)

izotopes = list(
  "углерод"= list(`12`= 0.9894, `13`= 0.0106),
  "водород"= list(`1`= 0.999855, `2`= 0.000145),
  "азот"= list(`14`= 0.996205, `15`= 0.003795),
  "кислород"= list(`16`= 0.9976, `17`= 0.0004, `18`= 0.002),
  "сера"= list(
    `32`= 0.9485, 
    `33`= 0.00763, 
    `34`= 0.04365, 
    `35`= 0.0,      
    `36`= 0.000158
  ),
  "хлор"= list(
    `35`= 0.758, 
    `36`= 0.0,     
    `37`= 0.242
  ),
  "бром"= list(
    `79`= 0.5065, 
    `80`= 0.0,      
    `81`= 0.4935
  ),
  "фтор"=list(
    `19`=1
  ),
  "фосфор"=list(
    `31`=1
  )
)

exact_convolve <- function(x, y) {
  m <- outer(x, y)
  indices <- outer(seq_along(x), seq_along(y), "+")
  res <- as.vector(tapply(m, indices, sum))
  return(res)
}

# Exact molecular weight
mol_weight_calc <- function(atoms) {
  molecular_weight = 0
  for (i in names(atoms)) {
    if (atoms[[i]] > 0) {
      for (j in names(izotopes_1[[i]])) {
        molecular_weight = molecular_weight + atoms[[i]] * izotopes_1[[i]][[j]] * as.numeric(j)
      }
    }
  }
  return(molecular_weight)
}

# Molecular mass of the most abundant
mol_weight_abundant <- function(atoms) {
  mol_weight_abundant = 0
  for (i in names(atoms)) {
    if (atoms[[i]] > 0) {
      k = names(izotopes[[i]])
      mol_weight_abundant = mol_weight_abundant + atoms[[i]] * as.numeric(k[1])
    }
  }
  return(as.integer(mol_weight_abundant))
}

# Isotopic distribution through polinomial folding 

calc_isotope_distribution <- function(atoms) {
  total_dist <- c("0" = 1.0)
  
  for (element in names(atoms)) {
    count <- as.integer(atoms[[element]])
    if (is.na(count) || count <= 0) next
    
    elem_iso <- izotopes[[element]]
    iso_masses <- as.integer(names(elem_iso))
    iso_probs <- unlist(elem_iso, use.names = FALSE)
    
    min_em <- min(iso_masses)
    max_em <- max(iso_masses)
    full_em_masses <- min_em:max_em
    
    single_atom_dist <- numeric(length(full_em_masses))
    names(single_atom_dist) <- as.character(full_em_masses)
    single_atom_dist[as.character(iso_masses)] <- iso_probs
    
    elem_dist <- single_atom_dist
    
    if (count > 1) {
      for (i in 2:count) {
        m1 <- as.integer(names(elem_dist))
        m2 <- as.integer(names(single_atom_dist))
        
        new_masses <- (min(m1) + min(m2)):(max(m1) + max(m2))
        
        elem_dist <- exact_convolve(elem_dist, single_atom_dist)
        names(elem_dist) <- as.character(new_masses)
      }
    }
    
    m_total <- as.integer(names(total_dist))
    m_elem <- as.integer(names(elem_dist))
    
    comb_masses <- (min(m_total) + min(m_elem)):(max(m_total) + max(m_elem))
    
    total_dist <- exact_convolve(total_dist, elem_dist)
    names(total_dist) <- as.character(comb_masses)
  }
  
  return(total_dist)
}


find_percent_of_consist_fast <- function(atoms, mol_mass_needed) {
  dist <- calc_isotope_distribution(atoms)
  mass_str <- as.character(mol_mass_needed)
  
  if (mass_str %in% names(dist)) {
    return(unname(dist[mass_str]))
  } else {
    return(0.0)
  }
}

# ==== UI ====
ui <- page_sidebar(
  theme = bs_theme(
    version = 5, 
    bootswatch = "lux",                
    base_font = font_google('Inter')    
  ) %>%
    bs_add_rules(
      "
      /* CSS */

body {
  background-color: #fcfcfd;
  font-family: -apple-system, BlinkMacSystemFont, 'Inter', sans-serif;
}

.sidebar {
  background-color: #ffffff !important;
  border-right: 1px solid #f1f1f4;
  box-shadow: none;
}

.card {
  border: 1px solid #eef0f3;
  background-color: #ffffff;
  box-shadow: 0 1px 3px rgba(0,0,0,0.02);
  border-radius: 10px;
  transition: all 0.25s cubic-bezier(0.4, 0, 0.2, 1);
}
.card:hover {
  transform: translateY(-4px);
  box-shadow: 0 12px 24px rgba(0,0,0,0.04);
  border-color: #dcdfe4;
}

.result-title {
  font-size: 0.8rem;
  text-transform: uppercase;
  letter-spacing: 1px;
  color: #8a92a6;
  font-weight: 600;
}
.result-value {
  font-size: 1.8rem;
  font-weight: 600;
  color: #111827;
}

.formula-badge {
  font-size: 1.1rem;
  padding: 8px 16px;
  background-color: #f4f4f7;
  color: #111827;
  border-radius: 6px;
  font-weight: 500;
  border: 1px solid #e4e4e7;
}

.accent-card {
  background: #111827 !important;
  color: white !important;
  border: none;
}
.accent-card .result-title { color: #9ca3af; }
.accent-card .result-value { color: #ffffff; font-size: 2.2rem; }

    "),
  
  title = "Isotope Ratio & Mass Calculator",
  
  sidebar = sidebar(
    title = "Elements & Parameters",
    width = 320,
    numericInput('carbon', 'Number of Carbon (C)', value = 5, min = 0, max = 1000),
    numericInput('hydrogen', 'Number of Hydrogen (H)', value = 10, min = 0, max = 1000),
    numericInput('nitrogen', 'Number of Nitrogen (N)', value = 0, min = 0, max = 1000),
    numericInput('oxygen', 'Number of Oxygen (O)', value = 0, min = 0, max = 1000),
    numericInput('sulfur', 'Number of Sulfur (S)', value = 0, min = 0, max = 1000),
    numericInput('chlorine', 'Number of Chlorine (Cl)', value = 0, min = 0, max = 1000),
    numericInput('bromine', 'Number of Bromine (Br)', value = 0, min = 0, max = 1000),
    numericInput('fluorine', 'Number of Fluorine (F)', value = 0, min = 0, max = 1000),
    numericInput('phosphorus', 'Number of Phosphorus (P)', value = 0, min = 0, max = 1000),
  ),
  
  layout_column_wrap(
    width = 1,
    layout_column_wrap(width = 1 / 2,
                       # Mass difference card
                       card(card_body(
                         div(class = "result-title", "Mass Difference (Δm)"),
                         div(style = "margin-top: 10px;", numericInput(
                           'diff',
                           '',
                           value = 3,
                           min = 0,
                           max = 1000
                         ))
                       )),
                       # Formula card
                       card(card_body(
                         div(class = "result-title", "Current Molecular Formula"),
                         div(style = "margin-top: 10px;", uiOutput("formula_ui"))
                       ))),
    
    layout_column_wrap(
      width = 1 / 3,
      
      card(card_body(
        div(class = "result-title", "Average Molecular Mass"),
        div(class = "result-value", textOutput("mol_mass", inline = TRUE))
      )),
      
      card(card_body(
        div(class = "result-title", "Monoisotopic / Needed Mass"),
        div(class = "result-value", textOutput("mol_mas_abundant", inline = TRUE))
      )),
      card(class = "accent-card",
           card_body(
             div(class = "result-title", "Abundance Percentage"),
             div(class = "result-value", textOutput("percent", inline = TRUE))
           ))
    )
  )
)

# ==== Server ====
server <- function(input, output) {
  
  atoms <- reactive({
    list(
      "углерод" = input$carbon,
      "водород" = input$hydrogen,
      "азот"    = input$nitrogen,
      "кислород"= input$oxygen,
      "сера"    = input$sulfur,
      "хлор"    = input$chlorine,
      "бром"    = input$bromine,
      "фтор"    = input$fluorine,
      "фосфор"  = input$phosphorus
    )
  })
  
  molacular_mass_needed <- reactive({
    mol_weight_abundant(atoms()) + as.integer(input$diff)
  })
  
  # Formula
  output$formula_ui <- renderUI({
    parts <- list()
    elements <- list("C" = input$carbon, "H" = input$hydrogen, "N" = input$nitrogen, 
                     "O" = input$oxygen, "S" = input$sulfur, "Cl" = input$chlorine, "Br" = input$bromine,
                     "F" = input$fluorine, "P" = input$phosphorus)
    
    for (el in names(elements)) {
      val <- as.integer(elements[[el]])
      if (!is.na(val) && val > 0) {
        parts <- c(parts, list(el, tags$sub(val)))
      }
    }
    
    if (length(parts) == 0) return(span(class = "formula-badge", "Empty"))
    do.call(div, c(list(class = "formula-badge"), parts))
  })
  
  output$mol_mass <- renderText({
    req(atoms())
    sprintf("%.5f Da", mol_weight_calc(atoms()))
  })
  
  output$mol_mas_abundant <- renderText({
    req(atoms())
    paste0(mol_weight_abundant(atoms()), " → ", molacular_mass_needed(), " Da")
  })
  
  output$percent <- renderText({
    req(atoms(), molacular_mass_needed())
    res_percent <- find_percent_of_consist_fast(atoms(), molacular_mass_needed())
    sprintf("%.6f %%", res_percent * 100)
  })
}

shinyApp(ui = ui, server = server)
