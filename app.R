library(shiny)
library(purrr)
library(dplyr)
library(DT)
library(shinythemes)
library(ggplot2)
library(iNEXT)
library(vegan)
library(hillR)

# Load ggplot theme

theme_thomas <- function(){
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'grey30'),
        
        panel.grid.major = element_line(color = "grey92", 
                                        linewidth = .4),
        
        panel.grid.minor = element_blank(),
        
        axis.title.x = element_text(color = "grey30", 
                                    margin = margin(t = 20), 
                                    size = 24),
        
        axis.title.y = element_text(color = "grey30", 
                                    margin = margin(r = 20), 
                                    size = 24),
        
        axis.text = element_text(color = "grey50", 
                                 size = 20),
        
        axis.ticks =  element_line(color = "grey92", 
                                   linewidth = .5),
        
        axis.ticks.length = unit(.6, "lines"),
        
        legend.position = "top",
        
        plot.title = element_text(hjust = 0, 
                                  color = "black", 
                                  family = "Roboto",
                                  size = 21, 
                                  margin = margin(t = 10, b = 35)),
        
        plot.subtitle = element_text(hjust = 0, 
                                     color = "grey30",
                                     family = "Roboto", 
                                     size = 14, 
                                     margin = margin(0, 0, 25, 0)),
        
        plot.title.position = "plot",
        
        plot.caption = element_text(color = "grey50", 
                                    size = 10, 
                                    hjust = 1,
                                    family = "Roboto", 
                                    lineheight = 1.05, 
                                    margin = margin(30, 0, 0, 0)),
        
        plot.caption.position = "plot", 
        
        plot.margin = margin(rep(20, 4)), 
        
        panel.spacing=unit(4,"lines"))
  
}

# Create the data:

# Set seed for reproducibility
set.seed(123)

# Species names
species_names <- c("Adelpha attica", "Adelpha boeotia", "Adelpha capucinus", "Adelpha cocala", "Adelpha delinita",
                   "Adelpha erotia", "Adelpha heraclea", "Adelpha iphiclus", "Adelpha jordani", "Adelpha malea",
                   "Adelpha mesentina", "Adelpha messana", "Adelpha naxia", "Adelpha plesaure", "Adelpha pollina",
                   "Adelpha thesprotia", "Archaeoprepona demophon", "Archaeoprepona demophoon", "Archaeoprepona licomedes",
                   "Baeotus aeilus", "Baeotus deucalion", "Bia actorion", "Caeruleuptychia scopulata", "Caligo eurilochus",
                   "Callicore hystaspes", "Catoblepia berecynthia", "Catoblepia soranus", "Catoblepia xanthicles",
                   "Catonephele acontius", "Catonephele numilia", "Chloreuptychia callichloris", "Chloreuptychia hersis",
                   "Cissia myneca", "Cissia terrestris", "Colobura annulata", "Colobura dirce", "Consul fabius",
                   "Epiphile lampethusa", "Erichthodes antonina", "Eunica caelina", "Eunica celma", "Eunica eurota",
                   "Eunica marsolia", "Eunica orphise", "Eunica sophonisba", "Eunica sydonia", "Eunica veronica",
                   "Euptychia ordinata", "Euptychia sp", "Fountainea ryphea ryphea", "Hamadryas arinome", "Hamadryas feronia",
                   "Hamadryas laodamia", "Harjesia n sp", "Harjesia obscura", "Historis acheronta", "Historis odius",
                   "Hypna clytemnestra", "Magneuptychia fugitiva", "Manataria hercyna", "Megeuptychia antonoe",
                   "Memphis acidalia", "Memphis basilia", "Memphis catinka", "Memphis glauce", "Memphis moruus",
                   "Memphis offa", "Memphis phantes vicinia", "Memphis philumena", "Memphis polycarmes", "Memphis xenocles",
                   "Memphis xxx", "Morpho achilles", "Morpho helenor", "Myscelia capenas", "Nessaea hewitsonii",
                   "Nessaea obrinus", "Opoptera aorsa", "Opsiphanes quiteria", "Panacea prola", "Panacea regina",
                   "Pareuptychia", "Paulogramma peristera", "Posttaygetis penelea", "Prepona laertes", "Pseudodebis celia",
                   "Pseudodebis valentina", "Pyrrhogyra amphiro", "Pyrrhogyra crameri", "Pyrrhogyra edocla", "Pyrrhogyra otolais",
                   "Smyrna blomfildia", "Splendeuptychia n sp", "Taygetis cleopatra/sosis", "Taygetis elegia",
                   "Taygetis laches/thamyra", "Taygetis n sp", "Taygetis oyapock", "Taygetis sylvia", "Taygetis virgilia",
                   "Taygetis xx", "Temenis laothoe", "Temenis pulchra", "Tigridia acesta", "Zaretis isidora",
                   "Zaretis itys", "Zaretis n sp", "Zaretis strigosus")


# Number of species for each dataset
num_species_old_growth <- 89
num_species_second_growth <- 54

# Ensure both datasets together have 108 species

shared_species_nr <- sum(num_species_old_growth, num_species_second_growth) - length(species_names)

shared_species <- sample(species_names, shared_species_nr, replace = FALSE)

old_growth_species_unique <- sample(x = setdiff(species_names, shared_species), 
                                    size = (num_species_old_growth - shared_species_nr), 
                                    replace = FALSE)

second_growth_species_unique <- setdiff(species_names, c(shared_species, old_growth_species_unique))

old_growth_community <- c(shared_species, old_growth_species_unique)
second_growth_community <- c(shared_species, second_growth_species_unique)


# Set the community relative abundance distribution

# Number of individuals for each dataset
num_individuals_old_growth <- 5000
num_individuals_second_growth <- 10000

# Number of common species for each dataset
num_common_species_old_growth <- 11
num_common_species_second_growth <- 26

# Create abundance vectors for old-growth and second-growth communities
abundance_old_growth <- numeric(num_species_old_growth)
abundance_second_growth <- numeric(num_species_second_growth)

# Assign abundance for common species in old-growth
common_species_old_growth <- sample(old_growth_community, num_common_species_old_growth, replace = FALSE)
abundance_old_growth[old_growth_community %in% common_species_old_growth] <- round(runif(num_common_species_old_growth, min = 500, max = 2000))

# Assign abundance for rare species in old-growth
abundance_old_growth[abundance_old_growth == 0] <- round(runif(sum(abundance_old_growth == 0), min = 1, max = 200))

# Assign abundance for common species in second-growth
common_species_second_growth <- sample(second_growth_community, num_common_species_second_growth, replace = FALSE)
abundance_second_growth[second_growth_community %in% common_species_second_growth] <- round(runif(num_common_species_second_growth, min = 500, max = 2000))

abundance_second_growth[abundance_second_growth == 0] <- round(runif(sum(abundance_second_growth == 0), min = 1, max = 500))

# Normalize abundance to achieve the desired total number of individuals
abundance_old_growth <- round(abundance_old_growth / sum(abundance_old_growth) * num_individuals_old_growth, digits = 0)
abundance_second_growth <- round(abundance_second_growth / sum(abundance_second_growth) * num_individuals_second_growth, digits = 0)

#Create datasets 

dataset_old_growth <- data.frame(species = old_growth_community, abundance = abundance_old_growth)
dataset_second_growth <- data.frame(species = second_growth_community, abundance = abundance_second_growth)

dataset_old_growth <- unlist(apply(dataset_old_growth, MARGIN = 1, function(x) paste0(rep(x[1], x[2]),"_", seq(1, x[2], 1))))
dataset_old_growth <- sample(dataset_old_growth, length(dataset_old_growth), replace = FALSE)

dataset_second_growth <- unlist(apply(dataset_second_growth, MARGIN = 1, function(x) paste0(rep(x[1], x[2]),"_", seq(1, x[2], 1))))
dataset_second_growth <- sample(dataset_second_growth, length(dataset_second_growth), replace = FALSE)

# Define UI
ui <- fluidPage(
    theme = shinytheme("united"),
    titlePanel("NATF320: Biodiversity sampling practical"),
    navbarPage(
        "Butterfly",
        tabPanel("Introduction",
                 h3("Assignment Introduction"),
                 p("This is the introduction text for the assignment. You can add more information here."),
                 # Add more content as needed
        ),
        tabPanel("Sample",
                 sidebarLayout(
                     sidebarPanel(
                         selectInput("forest_type", "Select Forest Type", c("Old-Growth", "Second-Growth", "Both")),
                         actionButton("sample_button", "Draw Random Sample"),
                         actionButton("reset_button", "Reset"),
                         br(), 
                         br(), 
                         textOutput("counter_text"), 
                         textOutput("sampled_text_1"), 
                         textOutput("sampled_text_2"), 
                         textOutput("sampled_text_3")
                         
                     ),
                     mainPanel(plotOutput("butterfly_plot_1"),
                               br(),
                               textOutput("oldgrowth_header"),
                               br(),
                               DTOutput("sample_table_1"), 
                               br(),
                               textOutput("secondary_header"),
                               br(),
                               DTOutput("sample_table_2")
                               )
                 )
        ),
         tabPanel("Explore",
            tabsetPanel(
                tabPanel("Diversity Metrics",
                    h3("Explore Diversity Metrics"),
                    
                    br(),
                    
                    h4(strong("A. Traditional diversity metrics")),
                    br(),
                    textOutput("oldgrowth_header_2"),
                    br(),
                    DTOutput("traditional_table"), 
                    br(), 
                    textOutput("secondary_header_2"),
                    br(),
                    DTOutput("traditional_table_2"),
                    br(),
                    h4(strong("B. Hill-based diversity metrics")),
                    br(), 
                    textOutput("oldgrowth_header_3"),
                    br(),
                    DTOutput("hill_table"),
                    br(), 
                    textOutput("secondary_header_3"),
                    br(),
                    DTOutput("hill_table_2"),
                    br()
                ),
                tabPanel("Sampling procedure",
                         h2("Explore sampling procedure"),
                         br(),
                         h3("1. What do diversity values look like per sample?"),
                         br(),
                         h4(strong(em("1A. Number of individuals per sample"))),
                         br(),
                         plotOutput("individual_sample_plot"),
                         br(), 
                         h4(strong(em("1B. Number of species per sample"))),
                         br(), 
                         plotOutput("species_sample_plot"), 
                         br(),
                         h3("2. How do the cumulative diversity values change with an increasing number of samples?"),
                         br(), 
                         h4(strong(em("2A. Cumulative species richness with an increasing number of samples"))),
                         br(), 
                         plotOutput("species_by_sample_plot"), 
                         br(),
                         h4(strong(em("2B. Cumulative Shannon diversity with an increasing number of samples"))),
                         br(), 
                         plotOutput("shannon_by_samples_plot"), 
                         br(),
                         h4(strong(em("2C. Cumulative Simpson diversity with an increasing number of samples"))),
                         br(),
                         plotOutput("simpson_by_samples_plot"), 
                         br(),
                         h4(strong(em("2D. Cumulative evenness with an increasing number of samples"))),
                         br(),
                         plotOutput("evenness_by_samples_plot"), 
                         br(),
                         h4(strong(em("2E. Cumulative diversity (q = 1) with an increasing number of samples"))),
                         br(),
                         plotOutput("hill_q1_by_samples_plot"), 
                         br(), 
                         h4(strong(em("2F. Cumulative diversity (q = 2) with an increasing number of samples"))),
                         br(),
                         plotOutput("hill_q2_by_samples_plot"), 
                         br(), 
                         h4(strong(em("2G. Cumulative evenness (Hill-based) with an increasing number of samples"))),
                         br(),
                         plotOutput("hill_evenness_by_samples_plot"), 
                         br(),
                         h3("3. How do the cumulative diversity values change with an increasing number of sampled individuals?"),
                         br(),
                         h4(strong(em("3A. Cumulative species richness with an increasing number of individuals"))),
                         br(), 
                         plotOutput("species_by_individuals_plot"), 
                         br(),
                         h4(strong(em("3B. Cumulative Shannon diversity with an increasing number of individuals"))),
                         br(), 
                         plotOutput("shannon_by_individuals_plot"), 
                         br(),
                         h4(strong(em("3C. Cumulative Simpson diversity with an increasing number of individuals"))),
                         br(),
                         plotOutput("simpson_by_individuals_plot"), 
                         br(),
                         h4(strong(em("3D. Cumulative evenness with an increasing number of individuals"))),
                         br(),
                         plotOutput("evenness_by_individuals_plot"), 
                         br(),
                         h4(strong(em("3E. Cumulative diversity (q = 1) with an increasing number of individuals"))),
                         br(),
                         plotOutput("hill_q1_by_individuals_plot"), 
                         br(), 
                         h4(strong(em("3F. Cumulative diversity (q = 2) with an increasing number of individuals"))),
                         br(),
                         plotOutput("hill_q2_by_individuals_plot"), 
                         br(), 
                         h4(strong(em("3G. Cumulative evenness (Hill-based) with an increasing number of individuals"))),
                         br(),
                         plotOutput("hill_evenness_by_individuals_plot"), 
                         br()

                ),
                tabPanel("Abundance Distribution",
                    h3("Explore Abundance Distribution"),
                    plotOutput("butterfly_plot_2")
                    # ... (add UI elements specific to abundance distribution)
                ),
                tabPanel("Rarefaction plots",
                         h3("Explore rarefaction plots"),
                         br(), 
                         sidebarLayout(
                           sidebarPanel(
                             actionButton("rarefaction_button", "Rarefy data"),
                             actionButton("reset_button_rar", "Reset")), 
                           mainPanel(
                             h3("1. Size-based rarefaction"),
                             br(), 
                             h4(strong(em("1A. Size-based rarefaction plot (abundance)"))),
                             br(),
                             plotOutput("rar_plot_abun_siz"), 
                             br(), 
                             h4(strong(em("1B. Size-based rarefaction plot (incidence)"))),
                             br(),
                             plotOutput("rar_plot_inc_siz"), 
                             br(), 
                             h3("2. Coverage-based rarefaction"),
                             br(), 
                             h4(strong(em("2A. Coverage-based rarefaction plot (abundance)"))), 
                             br(),
                             plotOutput("rar_plot_abun_cov"), 
                             br(), 
                             h4(strong(em("2B. Coverage-based rarefaction plot (incidence)"))), 
                             plotOutput("rar_plot_inc_cov"))
                         )
                           ), 
                tabPanel("Rarefaction data", 
                         h3("Explore rarefaction results"),
                         br(),
                         h3("1. Size-based rarefaction"),
                         br(), 
                         h4(strong(em("1A. Size-based rarefaction table (abundance)"))),
                         br(),
                         DTOutput("rar_table_abun_siz"),
                         br(), 
                         h4(strong(em("1B. Size-based rarefaction table (incidence)"))),
                         br(),
                         DTOutput("rar_table_inc_siz"),
                         br(), 
                         h3("2. Coverage-based rarefaction"),
                         br(), 
                         h4(strong(em("2A. Coverage-based rarefaction table (abundance)"))), 
                         br(),
                         DTOutput("rar_table_abun_cov"),
                         br(), 
                         h4(strong(em("2B. Coverage-based rarefaction table (incidence)"))), 
                         DTOutput("rar_table_inc_cov")
                )
                         
            ),
        ),
        tabPanel("Compare",
                 # Add content for the Compare tab
                 h3("Compare Forest Types"),
                 # ... (add additional UI elements for Compare tab)
        ),
        tabPanel("Download",
                 # Add content for the Download tab
                 h3("Download Biodiversity Data"),
                 # ... (add additional UI elements for Download tab)
        )
    )
)

# Define server
server <- function(input, output, session) {
    
    # Define datasets
    dataset_old_growth <- dataset_old_growth  # Insert the code for dataset_old_growth creation here
    dataset_second_growth <- dataset_second_growth  # Insert the code for dataset_second_growth creation here
    dataset_both <- vector("list")
    dataset_both[[1]] <- dataset_old_growth
    dataset_both[[2]] <- dataset_second_growth
    
    # 1. Set counter for the number of sample button clicks
    sample_counter <- reactiveVal(0)
    
    observeEvent(input$sample_button, {
        # Increment the counter
        sample_counter(sample_counter() + 1)
        
        # Update the UI element displaying the counter
        output$counter_text <- renderText({
            paste("Number of samples:", sample_counter())
        })

    })
    

    # 2. Make a reset button
    
    observeEvent(input$reset_button, {
        shinyjs::enable("sample_button")  # Re-enable the sample button
        sample_counter(0)
        
        # Reset the forest type to NULL
        updateSelectInput(session, "forest_type", selected = NULL)
        
        butterfly_dataset <- reactiveVal(NULL)
        butterfly_dataset_old <- reactiveVal(NULL)
        butterfly_dataset_sec <- reactiveVal(NULL)
        
        # Update the UI element displaying the sampled_text
        output$sampled_text_1 <- renderText({
            paste("Number of individuals sampled:", 0)
        })
        
        output$sampled_text_2 <- renderText({
          NULL
        })
        
        output$sampled_text_3 <- renderText({
          NULL
        })
        
        output$oldgrowth_header <- renderText({
          NULL
        })
        
        output$oldgrowth_header_2 <- renderText({
          NULL
        })
        
        output$oldgrowth_header_3 <- renderText({
          NULL
        })
        
        output$secondary_header <- renderText({
          NULL
        })
        
        output$secondary_header_2 <- renderText({
          NULL
        })
        
        output$secondary_header_3 <- renderText({
          NULL
        })

        
        # Clear the plot, table, and butterfly_data
        output$butterfly_plot_1 <- renderPlot(NULL)
        output$butterfly_plot_2 <- renderPlot(NULL)
        
        output$sample_table_1 <- renderDT(NULL)
        output$sample_table_2 <- renderDT(NULL)
        
        output$traditional_table <- renderDT(NULL)
        output$traditional_table_2 <- renderDT(NULL)
        
        output$hill_table <- renderDT(NULL)
        output$hill_table_2 <- renderDT(NULL)
        
        
        butterfly_data$data <- vector("list")
        butterfly_data_old$data <- vector("list")
        butterfly_data_sec$data <- vector("list")
        output$traditional_table <- renderDT(NULL)
        output$hill_table <- renderDT(NULL)
        
        output$rar_plot_abun_siz <- renderPlot(NULL)
        output$rar_plot_inc_siz <- renderPlot(NULL)
        output$rar_plot_abun_cov <- renderPlot(NULL)
        output$rar_plot_inc_cov <- renderPlot(NULL)
    })
    
    observeEvent(input$reset_button_rar, {
      
      output$rar_plot_abun_siz <- renderPlot(NULL)
      output$rar_plot_inc_siz <- renderPlot(NULL)
      output$rar_plot_abun_cov <- renderPlot(NULL)
      output$rar_plot_inc_cov <- renderPlot(NULL)
      
    })
    
    
    # 2. Observe whether the drop-down changes, if it does, reset everything
    
    observe({
        old_forest_type <- isolate(input$forest_type)
        
        if (!is.null(old_forest_type)) {
            shinyjs::enable("sample_button")  # Re-enable the sample button
            sample_counter(0)
            
            # Update the selected dataset based on the new forest type
            selected_dataset(
                switch(input$forest_type,
                       "Old-Growth" = dataset_old_growth,
                       "Second-Growth" = dataset_second_growth, 
                       "Both" = dataset_both)
            )
            
            butterfly_dataset <- reactiveVal(NULL)
            butterfly_dataset_old <- reactiveVal(NULL)
            butterfly_dataset_sec <- reactiveVal(NULL)
            
            # Update the UI element displaying the sampled_text
            output$sampled_text_1 <- renderText({
                paste("Number of individuals sampled:", 0)
            })
            
            output$sampled_text_2 <- renderText({
              NULL
            })
            
            output$sampled_text_3 <- renderText({
              NULL
            })
            
            output$oldgrowth_header <- renderText({
              NULL
            })
            
            output$oldgrowth_header_2 <- renderText({
              NULL
            })
            
            output$oldgrowth_header_3 <- renderText({
              NULL
            })
            
            
            
            output$secondary_header <- renderText({
              NULL
            })
            
            output$secondary_header_2 <- renderText({
              NULL
            })
            
            output$secondary_header_3 <- renderText({
              NULL
            })
            
            # Clear the plot and table
            output$butterfly_plot_1 <- renderPlot(NULL)
            output$butterfly_plot_2 <- renderPlot(NULL)
            
            output$sample_table_1 <- renderDT(NULL)
            output$sample_table_2 <- renderDT(NULL)
            
            output$traditional_table <- renderDT(NULL)
            output$traditional_table_2 <- renderDT(NULL)
            
            output$hill_table <- renderDT(NULL)
            output$hill_table_2 <- renderDT(NULL)
            
            butterfly_data$data <- vector("list")
            butterfly_data_old$data <- vector("list")
            butterfly_data_sec$data <- vector("list")
            
            output$traditional_table <- renderDT(NULL)
            output$hill_table <- renderDT(NULL)
            
            output$rar_plot_abun_siz <- renderPlot(NULL)
            output$rar_plot_inc_siz <- renderPlot(NULL)
            output$rar_plot_abun_cov <- renderPlot(NULL)
            output$rar_plot_inc_cov <- renderPlot(NULL)
        }
    }, priority = 1000)
    
    # 3. Create reactive object to save sampling data
    
    butterfly_data <- reactiveValues(data = list())
    butterfly_data_old <- reactiveValues(data = list())
    butterfly_data_sec <- reactiveValues(data = list())
    
    # 4. Create reactive object to store community data after sampling
    
    selected_dataset <- reactiveVal(NULL)
    num_individuals_sample <- reactiveVal(NULL)
    color_var <- reactiveVal(NULL)
    
    butterfly_dataset <- reactiveVal(NULL)
    butterfly_dataset_old <- reactiveVal(NULL)
    butterfly_dataset_sec <- reactiveVal(NULL)
    
    
    observeEvent(input$forest_type, {
        selected_dataset(
            switch(input$forest_type,
                   "Old-Growth" = dataset_old_growth,
                   "Second-Growth" = dataset_second_growth, 
                   "Both" = dataset_both)
        )
    })
    
    # 5. Create events when sample button is clicked 
    
    
    observeEvent(input$sample_button, {
        
        
        # 5.2. Sample the community
        
        if(input$forest_type == "Old-Growth" | input$forest_type == "Second-Growth"){
          
          if(input$forest_type == "Old-Growth"){
            
            num_individuals_sample(sample(c(1, 1, 2, 2, 2, 2, 3, 3), 
                                             size = 1))
            
            color_var("#769a6e")
            
          }
          
          else{
            
            if(input$forest_type == "Second-Growth"){
              
              num_individuals_sample(sample(c(1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 8, 8, 9, 9, 10, 11, 12), 
                                                size = 1))
              
              color_var("#FFBF00")
              
            }
        
          }
        
        sample_result <- sample(selected_dataset(), size = num_individuals_sample(), replace = FALSE)
        
        # 5.3. Update the selected dataset, removing the sampled individuals 
        
        selected_dataset(setdiff(selected_dataset(), sample_result))
        
        # 5.4. Prepare the sampling data
        
        sample_result_df <- data.frame(table(gsub("_.*","",sample_result)))
        colnames(sample_result_df) <- c("species", "sample")
        
        if (length(butterfly_data$data) == 0) {
            butterfly_data$data[[1]] <- sample_result_df
        } else {
            butterfly_data$data[[length(butterfly_data$data) + 1]] <- sample_result_df
        }
        
        # 5.5. Create output sample table from data

        butterfly_summary <- reduce(butterfly_data$data, full_join, by = "species")
        butterfly_summary[is.na(butterfly_summary)] <- 0
        rownames <- butterfly_summary[,1]
        butterfly_summary <- as.data.frame(butterfly_summary[,(2:ncol(butterfly_summary))])
        rownames(butterfly_summary) <- rownames
        colnames(butterfly_summary) <- paste0("S_", seq(1, length(butterfly_data$data), 1))
        
        butterfly_dataset(butterfly_summary)
        
        # 5.6. Create counter of how many individuals sampled
        
        output$sampled_text_1 <- renderText({
            paste("Number of individuals:", sum(unlist(butterfly_summary)))
        })
        
        output$sampled_text_2 <- renderText({
          NULL
        })
        
        output$sampled_text_3 <- renderText({
          NULL
        })
        
        output$oldgrowth_header <- renderText({
          NULL
        })
        
        output$oldgrowth_header_2 <- renderText({
          NULL
        })
        
        output$oldgrowth_header_3 <- renderText({
          NULL
        })
        
        output$secondary_header <- renderText({
          NULL
        })
        
        output$secondary_header_2 <- renderText({
          NULL
        })
        
        output$secondary_header_3 <- renderText({
          NULL
        })
        
        
        # 5.6. Display sample data table
        output$sample_table_1 <- renderDT({
            datatable(
                butterfly_summary,
                rownames = TRUE,
                options = list(scrollX = TRUE, scrollY = "400px")
            )
        })
        
        output$sample_table_2 <- renderDT({
          
          NULL
          
        })
        
        output$traditional_table_2 <- renderDT({
          
          NULL
          
        })
        
        output$hill_table_2 <- renderDT({
          
          NULL
          
        })
        
        # 5.7. Create dataset for rank abundance plot 
        
        rank_abundance_data <- data.frame(count = rowSums(butterfly_summary), 
                                   species = rownames(butterfly_summary))
        
        rank_abundance_data <- rank_abundance_data %>%
            mutate(rank = rank(desc(count), ties.method = "random"))
        
        # 5.8. Render rank abundance plot
        
        output$butterfly_plot_1 <- output$butterfly_plot_2 <- 
            
            renderPlot({
            
            ggplot(rank_abundance_data,
                   aes(x = rank,
                       y = count, 
                       label = species)) + 
                geom_line(color = color_var(),  
                          linewidth = 1) + 
                geom_point(size = 4, 
                           shape = 21, 
                           color = color_var(), 
                           fill = "white", 
                           stroke = 1) + 
                
                theme_bw()+
                
                geom_text(hjust=-0.1, vjust=0.1, angle = 45)+
                
                xlab("Species rank")+
                ylab("Abundance")+
                ggtitle("Rank abundance plot")+
                
                scale_x_continuous(expand = c(0.05, 0.5))+
                scale_y_continuous(expand = c(0.1, 1), 
                                   limits = c(0, (max(rank_abundance_data$count)+(max(rank_abundance_data$count)/7))))
            
        })
        
        # 5.9 Format data for species-by-sample plot
        
        incidence_plot <- butterfly_summary
        
        incidence_plot[incidence_plot > 0] <- 1
        
        sample_plot_data <- data.frame(count = colSums(incidence_plot), 
                                       samples = seq(1, length(colSums(butterfly_summary)), 1))
      
        
        # 5.10. Render species-by-sample plot
        
        output$species_sample_plot  <- 
          
          renderPlot({
            
            ggplot(sample_plot_data,
                   aes(x = samples,
                       y = count)) + 
              geom_line(color = "black",  
                        linewidth = 1) + 
              geom_point(size = 4, 
                         shape = 21, 
                         color = "black", 
                         fill = "white", 
                         stroke = 1) + 
              theme_thomas()+ 
              
              xlab("Sample number")+
              ylab("Species richness")+
              
              scale_x_continuous(expand = c(0.05, 0.5))+
              scale_y_continuous(expand = c(0.1, 1))
            
          })
        
        
        # 5.9 Format data for individuals-by-sample plot
        
        sample_individual_data <- data.frame(count = colSums(butterfly_summary), 
                                       samples = seq(1, length(colSums(butterfly_summary)), 1))
        
        # 5.10. Render individual-by-sample plot
        
        output$individual_sample_plot  <- 
          
          renderPlot({
            
            ggplot(sample_individual_data,
                   aes(x = samples,
                       y = count)) + 
              geom_line(color = "black",  
                        linewidth = 1) + 
              geom_point(size = 4, 
                         shape = 21, 
                         color = "black", 
                         fill = "white", 
                         stroke = 1) + 
              theme_thomas()+ 
              
              xlab("Sample number")+
              ylab("Individuals")+
              
              scale_x_continuous(expand = c(0.05, 0.5))+
              scale_y_continuous(expand = c(0.1, 1))
            
          })
        
        # Format data for cumulative number of individuals
        
        individual_acc_list <- vector("list")
        
        if(ncol(butterfly_summary) > 1){
          
          for (i in 1:ncol(butterfly_summary)){
            
            individual_acc_list[[i]] <- butterfly_summary[,c(1:i)]
            
          } 
          
          for (i in 2:length(individual_acc_list)){
            
            individual_acc_list[[i]] <- rowSums(individual_acc_list[[i]])
            
          }
          
        }
        
        # Format data for cumulative number of species
        
        species_acc_list <- vector("list")
        
        if(ncol(butterfly_summary) > 1){
          
          for (i in 1:ncol(butterfly_summary)){
            
            species_acc_list[[i]] <- butterfly_summary[,c(1:i)]
            species_acc_list[[1]] <- sum(species_acc_list[[1]]>0)
            
          } 
          
          for (i in 2:length(species_acc_list)){
            
            species_acc_list[[i]] <- sum(rowSums(species_acc_list[[i]])>0)
            
          }
          
        }
        
        # 5.12. Render species/diversity accumulation plots 
        
        if(ncol(butterfly_summary) > 1){
          
          # Species by samples
          
          species_by_sample_df <- data.frame(sample = seq(1, ncol(butterfly_summary), 1), 
                                       count = unlist(species_acc_list))
          
          output$species_by_sample_plot  <- 
            
            renderPlot({
              
              ggplot(species_by_sample_df,
                     aes(x = sample,
                         y = count)) + 
                geom_line(color = "black",  
                          linewidth = 1) + 
                geom_point(size = 4, 
                           shape = 21, 
                           color = "black", 
                           fill = "white", 
                           stroke = 1) + 
                theme_thomas()+
                
                xlab("Sample number")+
                ylab("Species richness")+
                
                scale_x_continuous(expand = c(0.05, 0.5))+
                scale_y_continuous(expand = c(0.1, 1))
              
            })
          
          
          # Species by individuals
          
          species_by_individuals_df <- data.frame(count = unlist(species_acc_list), 
                                          individuals = sapply(individual_acc_list, function(x)
                                            sum(unlist(x))))
          
          output$species_by_individuals_plot  <- 
            
            renderPlot({
              
              ggplot(species_by_individuals_df,
                     aes(x = individuals,
                         y = count)) + 
                geom_line(color = "black",  
                          linewidth = 1) + 
                geom_point(size = 4, 
                           shape = 21, 
                           color = "black", 
                           fill = "white", 
                           stroke = 1) + 
                theme_thomas()+
                
                xlab("Number of individuals")+
                ylab("Species richness")+
                
                scale_x_continuous(expand = c(0.05, 0.5))+
                scale_y_continuous(expand = c(0.1, 1))
              
            })
          
          # Shannon by samples
          
    shannon_by_samples_df <- data.frame(sample = seq(1, ncol(butterfly_summary), 1), 
          count = sapply(individual_acc_list, function(x)
            vegan::diversity(unlist(x), index = "shannon")))
    
    output$shannon_by_samples_plot  <- 
      
      renderPlot({
        
        ggplot(shannon_by_samples_df,
               aes(x = sample,
                   y = count)) + 
          geom_line(color = "black",  
                    linewidth = 1) + 
          geom_point(size = 4, 
                     shape = 21, 
                     color = "black", 
                     fill = "white", 
                     stroke = 1) + 
          theme_thomas()+
          
          xlab("Sample number")+
          ylab("Shannon's H'")+
          
          scale_x_continuous(expand = c(0.05, 0.5))+
          scale_y_continuous(expand = c(0.1, 1))
        
      })
    
    # Shannon by individuals
    
    shannon_by_individuals_df <- data.frame(individuals = sapply(individual_acc_list, function(x)
      sum(unlist(x))), 
                                        count = sapply(individual_acc_list, function(x)
                                          vegan::diversity(unlist(x), index = "shannon")))
    
    output$shannon_by_individuals_plot  <- 
      
      renderPlot({
        
        ggplot(shannon_by_individuals_df,
               aes(x = individuals,
                   y = count)) + 
          geom_line(color = "black",  
                    linewidth = 1) + 
          geom_point(size = 4, 
                     shape = 21, 
                     color = "black", 
                     fill = "white", 
                     stroke = 1) + 
          theme_thomas()+
          
          xlab("Number of individuals")+
          ylab("Shannon's H'")+
          
          scale_x_continuous(expand = c(0.05, 0.5))+
          scale_y_continuous(expand = c(0.1, 1))
        
      })
    
    # Simpson by samples
    
    simpson_by_samples_df <- data.frame(sample = seq(1, ncol(butterfly_summary), 1), 
                                 count = sapply(individual_acc_list, function(x)
                                   vegan::diversity(unlist(x), index = "simpson")))
    
    output$simpson_by_samples_plot  <- 
      
      renderPlot({
        
        ggplot(simpson_by_samples_df,
               aes(x = sample,
                   y = count)) + 
          geom_line(color = "black",  
                    linewidth = 1) + 
          geom_point(size = 4, 
                     shape = 21, 
                     color = "black", 
                     fill = "white", 
                     stroke = 1) + 
          theme_thomas()+
          
          xlab("Sample number")+
          ylab("Simpson's D")+
          
          scale_x_continuous(expand = c(0.05, 0.5))+
          scale_y_continuous(expand = c(0.1, 1))
        
      })
    
    # Simpson by individuals
    
    simpson_by_individuals_df <- data.frame(individuals = sapply(individual_acc_list, function(x)
      sum(unlist(x))), 
                                        count = sapply(individual_acc_list, function(x)
                                          vegan::diversity(unlist(x), index = "simpson")))
    
    output$simpson_by_individuals_plot  <- 
      
      renderPlot({
        
        ggplot(simpson_by_individuals_df,
               aes(x = individuals,
                   y = count)) + 
          geom_line(color = "black",  
                    linewidth = 1) + 
          geom_point(size = 4, 
                     shape = 21, 
                     color = "black", 
                     fill = "white", 
                     stroke = 1) + 
          theme_thomas()+
          
          xlab("Number of individuals")+
          ylab("Simpson's D")+
          
          scale_x_continuous(expand = c(0.05, 0.5))+
          scale_y_continuous(expand = c(0.1, 1))
        
      })
    
    
    # Evenness by samples
    
    evenness_by_samples_df <- data.frame(sample = seq(1, ncol(butterfly_summary), 1), 
                                         count = sapply(individual_acc_list, function(x)
                                           (vegan::diversity(x = unlist(x), index = "shannon")/log(sum(unlist(x>0))))))
    
    
    output$evenness_by_samples_plot  <- 
      
      renderPlot({
        
        ggplot(evenness_by_samples_df,
               aes(x = sample,
                   y = count)) + 
          geom_line(color = "black",  
                    linewidth = 1) + 
          geom_point(size = 4, 
                     shape = 21, 
                     color = "black", 
                     fill = "white", 
                     stroke = 1) + 
          theme_thomas()+
          
          xlab("Sample number")+
          ylab("Pielou's Evenness")+
          
          scale_x_continuous(expand = c(0.05, 0.5))+
          scale_y_continuous(expand = c(0.1, 1))
        
      })
    
    # Evenness by individuals
    
    evenness_by_individuals_df <- data.frame(individuals = sapply(individual_acc_list, function(x)
      sum(unlist(x))), 
      count = sapply(individual_acc_list, function(x)
        (vegan::diversity(x = unlist(x), index = "shannon")/log(sum(unlist(x>0))))))
    
    
    output$evenness_by_individuals_plot <- 
      
      renderPlot({
        
        ggplot(evenness_by_individuals_df,
               aes(x = individuals,
                   y = count)) + 
          geom_line(color = "black",  
                    linewidth = 1) + 
          geom_point(size = 4, 
                     shape = 21, 
                     color = "black", 
                     fill = "white", 
                     stroke = 1) + 
          theme_thomas()+
          
          xlab("Number of individuals")+
          ylab("Pielou's Evenness")+
          
          scale_x_continuous(expand = c(0.05, 0.5))+
          scale_y_continuous(expand = c(0.1, 1))
        
      })
    
    
    # Diversity (q = 0) by samples
    
    hill_q1_by_samples_df <- data.frame(sample = seq(1, ncol(butterfly_summary), 1), 
                                 count = sapply(individual_acc_list, function(x)
                                   hillR::hill_taxa(comm = unlist(x), q = 1, MARGIN = 2)))
    
    
    output$hill_q1_by_samples_plot  <- 
      
      renderPlot({
        
        ggplot(hill_q1_by_samples_df,
               aes(x = sample,
                   y = count)) + 
          geom_line(color = "black",  
                    linewidth = 1) + 
          geom_point(size = 4, 
                     shape = 21, 
                     color = "black", 
                     fill = "white", 
                     stroke = 1) + 
          theme_thomas()+
          
          xlab("Sample number")+
          ylab("Diversity (q = 1)")+
          
          scale_x_continuous(expand = c(0.05, 0.5))+
          scale_y_continuous(expand = c(0.1, 1))
        
      })
    
    # Diversity (q = 0) by individuals
    
    hill_q1_by_individuals_df <- data.frame(individuals = sapply(individual_acc_list, function(x)
      sum(unlist(x))), 
                                        count = sapply(individual_acc_list, function(x)
                                          hillR::hill_taxa(comm = unlist(x), q = 1, MARGIN = 2)))
    
    
    output$hill_q1_by_individuals_plot  <- 
      
      renderPlot({
        
        ggplot(hill_q1_by_individuals_df,
               aes(x = individuals,
                   y = count)) + 
          geom_line(color = "black",  
                    linewidth = 1) + 
          geom_point(size = 4, 
                     shape = 21, 
                     color = "black", 
                     fill = "white", 
                     stroke = 1) + 
          theme_thomas()+
          
          xlab("Number of individuals")+
          ylab("Diversity (q = 1)")+
          
          scale_x_continuous(expand = c(0.05, 0.5))+
          scale_y_continuous(expand = c(0.1, 1))
        
      })
    
    # Diversity (q = 2) by samples
    
    hill_q2_by_samples_df <- data.frame(sample = seq(1, ncol(butterfly_summary), 1), 
                                 count = sapply(individual_acc_list, function(x)
                                   hillR::hill_taxa(comm = unlist(x), q = 2, MARGIN = 2)))
    
    output$hill_q2_by_samples_plot  <- 
      
      renderPlot({
        
        ggplot(hill_q2_by_samples_df,
               aes(x = sample,
                   y = count)) + 
          geom_line(color = "black",  
                    linewidth = 1) + 
          geom_point(size = 4, 
                     shape = 21, 
                     color = "black", 
                     fill = "white", 
                     stroke = 1) + 
          theme_thomas()+
          
          xlab("Sample number")+
          ylab("Diversity (q = 2)")+
          
          scale_x_continuous(expand = c(0.05, 0.5))+
          scale_y_continuous(expand = c(0.1, 1))
        
      })
    
    # Diversity (q = 2) by individuals
    
    hill_q2_by_individuals_df <- data.frame(individuals = sapply(individual_acc_list, function(x)
      sum(unlist(x))), 
                                        count = sapply(individual_acc_list, function(x)
                                          hillR::hill_taxa(comm = unlist(x), q = 2, MARGIN = 2)))
    
    output$hill_q2_by_individuals_plot  <- 
      
      renderPlot({
        
        ggplot(hill_q2_by_individuals_df,
               aes(x = individuals,
                   y = count)) + 
          geom_line(color = "black",  
                    linewidth = 1) + 
          geom_point(size = 4, 
                     shape = 21, 
                     color = "black", 
                     fill = "white", 
                     stroke = 1) + 
          theme_thomas()+
          
          xlab("Number of individuals")+
          ylab("Diversity (q = 2)")+
          
          scale_x_continuous(expand = c(0.05, 0.5))+
          scale_y_continuous(expand = c(0.1, 1))
        
      })
    
    # Evenness (Hill-based) by samples
    
    hill_evenness_by_samples_df <- data.frame(sample = seq(1, ncol(butterfly_summary), 1), 
                                       count = sapply(individual_acc_list, function(x)
                                         ((hillR::hill_taxa(comm = unlist(x), q = 2, MARGIN = 2)-1)/
                                            (hillR::hill_taxa(comm = unlist(x), q = 0, MARGIN = 2)-1))))
    
    output$hill_evenness_by_samples_plot  <- 
      
      renderPlot({
        
        ggplot(hill_evenness_by_samples_df,
               aes(x = sample,
                   y = count)) + 
          geom_line(color = "black",  
                    linewidth = 1) + 
          geom_point(size = 4, 
                     shape = 21, 
                     color = "black", 
                     fill = "white", 
                     stroke = 1) + 
          theme_thomas()+
          
          xlab("Sample number")+
          ylab("Evenness (Hill-based)")+
          
          scale_x_continuous(expand = c(0.05, 0.5))+
          scale_y_continuous(expand = c(0.1, 1))
        
      })
    
    # Evenness (Hill-based) by individuals
    
    hill_evenness_by_individuals_df <- data.frame(individuals = sapply(individual_acc_list, function(x)
      sum(unlist(x))), 
                                              count = sapply(individual_acc_list, function(x)
                                                ((hillR::hill_taxa(comm = unlist(x), q = 2, MARGIN = 2)-1)/
                                                   (hillR::hill_taxa(comm = unlist(x), q = 0, MARGIN = 2)-1))))
    
    output$hill_evenness_by_individuals_plot  <- 
      
      renderPlot({
        
        ggplot(hill_evenness_by_individuals_df,
               aes(x = individuals,
                   y = count)) + 
          geom_line(color = "black",  
                    linewidth = 1) + 
          geom_point(size = 4, 
                     shape = 21, 
                     color = "black", 
                     fill = "white", 
                     stroke = 1) + 
          theme_thomas()+
          
          xlab("Number of individuals")+
          ylab("Evenness (Hill-based)")+
          
          scale_x_continuous(expand = c(0.05, 0.5))+
          scale_y_continuous(expand = c(0.1, 1))
        
      })
          
        }
        
        
        
        
        # 5.13. Calculate the traditional diversity metrics
        
        if (ncol(butterfly_summary) > 1){
            
            sample_abundance_list <- vector("list")
            
            for (i in 1:ncol(butterfly_summary)){
                
                sample_abundance_list[[i]] <- butterfly_summary[,c(1:i)]
                
            } 
            
            for (i in 2:length(species_acc_list)){
                
                sample_abundance_list[[i]] <- rowSums(sample_abundance_list[[i]])
                
            }
        
            # Species richness
        
        species_richness <- lapply(sample_abundance_list, 
                                   function(x) sum(x>0))
        
            # Shannon index
            
        shannon_index <- lapply(sample_abundance_list, 
                                function(x) round(vegan::diversity(x = x, index = "shannon"),digits = 3))

            
            # Simpson index
            
            simpson_index <- lapply(sample_abundance_list, 
                                    function(x) round(vegan::diversity(x = x, index = "simpson"), digits = 3))

        
            # Species evenness
        
            species_evenness <- lapply(sample_abundance_list, 
                                       function(x) round((vegan::diversity(x = x, index = "shannon")/log(sum(x>0))), digits = 3))

            # q = 0
            
            hill_q0 <- lapply(sample_abundance_list, 
                              function(x) round(hillR::hill_taxa(comm = x, 
                                                                  q = 0, 
                                                                  MARGIN = 2), digits = 3))
            
            # q = 1
            
            hill_q1 <- lapply(sample_abundance_list, 
                              function(x) round(hillR::hill_taxa(comm = x, 
                                                                  q = 1, 
                                                                  MARGIN = 2), digits = 3))
                    
            
            # q = 2
            
            hill_q2 <- lapply(sample_abundance_list, 
                              function(x) round(hillR::hill_taxa(comm = x,
                                                                         q = 2, 
                                                                         MARGIN = 2), digits = 3))
            
            # Hill-based evenness
            
            hill_evenness <- lapply(sample_abundance_list, 
                                    function(x) round(((hillR::hill_taxa(comm = x, 
                                                                 q = 2, 
                                                                 MARGIN = 2)-1)) / (hillR::hill_taxa(comm = x, 
                                                                                                     q = 0, 
                                                                                                     MARGIN = 2)-1), 
                                              digits = 3))
            
            
            traditional_table <- rbind(unlist(species_richness), 
                                       unlist(shannon_index), 
                                       unlist(simpson_index), 
                                       unlist(species_evenness))
            
            rownames(traditional_table) <- c("Species richness", 
                                             "Shannon's H'", 
                                             "Simpson's D", 
                                             "Pielou's evenness")
            
            colnames(traditional_table) <- paste0("Sample_", 
                                                  seq(1, ncol(butterfly_summary), 1))
            
            hill_table <- rbind(unlist(hill_q0), 
                                unlist(hill_q1), 
                                unlist(hill_q2), 
                                unlist(hill_evenness))
            
            rownames(hill_table) <- c("Diversity (q = 0)", 
                                      "Diversity (q = 1)", 
                                      "Diversity (q = 2)", 
                                      "Hill evenness")
            
            colnames(hill_table) <- paste0("Sample_", 
                                                  seq(1, ncol(butterfly_summary), 1))
            
            
            output$traditional_table <- renderDT({
                datatable(
                    traditional_table,
                    rownames = TRUE,
                    options = list(scrollX = TRUE, scrollY = "250px")
                )
            })
            
            
            output$hill_table <- renderDT({
                datatable(
                    hill_table,
                    rownames = TRUE,
                    options = list(scrollX = TRUE, scrollY = "250px")
                )
            })
            
            
            
        }
        
        } 
        
        
        else{
          
          if(input$forest_type == "Both"){
              
              num_individuals_sample_old <- sample(c(1, 1, 2, 2, 2, 2, 3, 3), 
                                                   size = 1)
              
              num_individuals_sample_sec <- sample(c(1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 8, 8, 9, 9, 10, 11, 12), 
                                                   size = 1)
            
            
            
            sample_result_old <- sample(selected_dataset()[[1]], size = num_individuals_sample_old, replace = FALSE)
            sample_result_sec <- sample(selected_dataset()[[2]], size = num_individuals_sample_sec, replace = FALSE)
            
            # 5.3. Update the selected dataset, removing the sampled individuals 
            
            newdata_old <- setdiff(selected_dataset()[[1]], sample_result_old)
            newdata_sec <- setdiff(selected_dataset()[[2]], sample_result_sec)
            newdata <- vector("list")
            newdata[[1]] <- newdata_old
            newdata[[2]] <- newdata_sec
            
            selected_dataset(newdata)


            # 5.4. Prepare the sampling data

            sample_result_df_old <- data.frame(table(gsub("_.*","",sample_result_old)))
            colnames(sample_result_df_old) <- c("species", "sample")
            
            sample_result_df_sec <- data.frame(table(gsub("_.*","",sample_result_sec)))
            colnames(sample_result_df_sec) <- c("species", "sample")

            if (length(butterfly_data_old$data) == 0) {
              butterfly_data_old$data[[1]] <- sample_result_df_old
            } else {
              butterfly_data_old$data[[length(butterfly_data_old$data) + 1]] <- sample_result_df_old
            }
            
            if (length(butterfly_data_sec$data) == 0) {
              butterfly_data_sec$data[[1]] <- sample_result_df_sec
            } else {
              butterfly_data_sec$data[[length(butterfly_data_sec$data) + 1]] <- sample_result_df_sec
            }

            # 5.5. Create output sample table from data

            butterfly_summary_old <- reduce(butterfly_data_old$data, full_join, by = "species")
            butterfly_summary_old[is.na(butterfly_summary_old)] <- 0
            rownames_old <- butterfly_summary_old[,1]
            butterfly_summary_old <- as.data.frame(butterfly_summary_old[,(2:ncol(butterfly_summary_old))])
            rownames(butterfly_summary_old) <- rownames_old
            colnames(butterfly_summary_old) <- paste0("S_", seq(1, length(butterfly_data_old$data), 1))
            
            butterfly_summary_sec <- reduce(butterfly_data_sec$data, full_join, by = "species")
            butterfly_summary_sec[is.na(butterfly_summary_sec)] <- 0
            rownames_sec <- butterfly_summary_sec[,1]
            butterfly_summary_sec <- as.data.frame(butterfly_summary_sec[,(2:ncol(butterfly_summary_sec))])
            rownames(butterfly_summary_sec) <- rownames_sec
            colnames(butterfly_summary_sec) <- paste0("S_", seq(1, length(butterfly_data_sec$data), 1))
            
            butterfly_dataset_old(butterfly_summary_old)
            butterfly_dataset_sec(butterfly_summary_sec)

            # 5.6. Create counter of how many individuals sampled
            
            output$sampled_text_1 <- renderText({
              
              NULL
              
            })

            output$sampled_text_2 <- renderText({
              paste("Number of individuals in old growth:", sum(unlist(butterfly_summary_old)))
            })
            
            output$sampled_text_3 <- renderText({
              paste("Number of individuals in secondary growth:", sum(unlist(butterfly_summary_sec)))
            })


            # 5.6. Display sample data table
            
            output$oldgrowth_header <- renderText({
              paste("Old-growth forest")
            })
            
            output$oldgrowth_header_2 <- renderText({
              paste("Old-growth forest")
            })
            
            output$oldgrowth_header_3 <- renderText({
              paste("Old-growth forest")
            })
            
            output$secondary_header <- renderText({
              paste("Secondary-growth forest")
            })
            
            output$secondary_header_2 <- renderText({
              paste("Secondary-growth forest")
            })
            
            output$secondary_header_3 <- renderText({
              paste("Secondary-growth forest")
            })
            
            output$sample_table_1 <-  renderDT({
              datatable(
                butterfly_summary_old,
                rownames = TRUE,
                options = list(scrollX = TRUE, scrollY = "400px"))
              
            })
            
            output$sample_table_2 <- renderDT({
              
              datatable(
              butterfly_summary_sec,
              rownames = TRUE,
              options = list(scrollX = TRUE, scrollY = "400px"))
              
            })
            
            
            # 5.7. Create dataset for rank abundance plot 
            
            rank_abundance_data_old <- data.frame(count = rowSums(butterfly_summary_old), 
                                              species = rownames(butterfly_summary_old))
            
            rank_abundance_data_old <- rank_abundance_data_old %>%
              mutate(rank = rank(desc(count), ties.method = "random"))
            
            rank_abundance_data_old$forest <- "Old"
            
            rank_abundance_data_sec <- data.frame(count = rowSums(butterfly_summary_sec), 
                                                  species = rownames(butterfly_summary_sec))
            
            rank_abundance_data_sec <- rank_abundance_data_sec %>%
              mutate(rank = rank(desc(count), ties.method = "random"))
            
            rank_abundance_data_sec$forest <- "Secondary"
            
            rank_abundance_data <- rbind(rank_abundance_data_old, 
                                         rank_abundance_data_sec)
            
            # 5.8. Render rank abundance plot
            
            output$butterfly_plot_1 <- output$butterfly_plot_2 <-

              renderPlot({

                ggplot(rank_abundance_data,
                       aes(x = rank,
                           y = count,
                           color = forest, 
                           fill = forest)) +
                  geom_line(linewidth = 1, 
                            alpha = 0.5) +
                  geom_point(size = 4,
                             shape = 21,
                             fill = "white",
                             stroke = 1, 
                             alpha = 0.5) +
                  
                  scale_color_manual(values = c("#769a6e", "#FFBF00"))+

                  theme_bw()+

                  xlab("Species rank")+
                  ylab("Abundance")+
                  ggtitle("Rank abundance plot")+

                  scale_x_continuous(expand = c(0.05, 0.5))+
                  scale_y_continuous(expand = c(0.1, 1),
                                     limits = c(0, (max(rank_abundance_data$count)+(max(rank_abundance_data$count)/7))))

              })
            
            
            # 5.9. Start generating biodiversity sampling plots
            
            incidence_plot_old <- butterfly_summary_old
            
            incidence_plot_old[incidence_plot_old > 0] <- 1
            
            sample_plot_data_old <- data.frame(count = colSums(incidence_plot_old), 
                                           samples = seq(1, length(colSums(butterfly_summary_old)), 1))
            
            sample_plot_data_old$forest_type <- "Old-growth"
            
            
            incidence_plot_sec <- butterfly_summary_sec
            
            incidence_plot_sec[incidence_plot_sec > 0] <- 1
            
            sample_plot_data_sec <- data.frame(count = colSums(incidence_plot_sec), 
                                               samples = seq(1, length(colSums(butterfly_summary_sec)), 1))
            
            sample_plot_data_sec$forest_type <- "Second-growth"
            
            sample_plot_data <- rbind(sample_plot_data_old, 
                                      sample_plot_data_sec)
            
            
            # 5.10. Render species-by-sample plot
            
            output$species_sample_plot  <-

              renderPlot({

                ggplot(sample_plot_data,
                       aes(x = samples,
                           y = count, 
                           color = forest_type)) +
                  geom_line(linewidth = 1) +
                  geom_point(size = 4,
                             shape = 21,
                             fill = "white",
                             stroke = 1) +
                  theme_thomas()+

                  xlab("Sample number")+
                  ylab("Species richness")+

                  scale_x_continuous(expand = c(0.05, 0.5))+
                  scale_y_continuous(expand = c(0.1, 1)) +
                  scale_color_manual(values = c("#769a6e", "#FFBF00"))

              })


            # 5.9 Format data for individuals-by-sample plot
            
            sample_individual_data_old <- data.frame(count = colSums(butterfly_summary_old),
                                                 samples = seq(1, length(colSums(butterfly_summary_old)), 1))
            
            sample_individual_data_old$forest_type <- "Old-growth"
            
            sample_individual_data_sec <- data.frame(count = colSums(butterfly_summary_sec),
                                                     samples = seq(1, length(colSums(butterfly_summary_sec)), 1))
            
            sample_individual_data_sec$forest_type <- "Second-growth"
            
            sample_individual_data <- rbind(sample_individual_data_old, 
                                            sample_individual_data_sec)
            
            # 5.10. Render individual-by-sample plot
            
            output$individual_sample_plot  <-

              renderPlot({

                ggplot(sample_individual_data,
                       aes(x = samples,
                           y = count, 
                           color = forest_type)) +
                  geom_line(linewidth = 1) +
                  geom_point(size = 4,
                             shape = 21,
                             fill = "white",
                             stroke = 1) +
                  theme_thomas()+

                  xlab("Sample number")+
                  ylab("Individuals")+

                  scale_x_continuous(expand = c(0.05, 0.5))+
                  scale_y_continuous(expand = c(0.1, 1)) + 
                  scale_color_manual(values = c("#769a6e", "#FFBF00"))

              })
            
            # Format data for cumulative number of individuals
            
            individual_acc_list_old <- vector("list")
            individual_acc_list_sec <- vector("list")
            
            if(ncol(butterfly_summary_old) > 1){
              
              for (i in 1:ncol(butterfly_summary_old)){
                
                individual_acc_list_old[[i]] <- butterfly_summary_old[,c(1:i)]
                individual_acc_list_sec[[i]] <- butterfly_summary_sec[,c(1:i)]
                
              } 
              
              for (i in 2:length(individual_acc_list_old)){
                
                individual_acc_list_old[[i]] <- rowSums(individual_acc_list_old[[i]])
                individual_acc_list_sec[[i]] <- rowSums(individual_acc_list_sec[[i]])
                
              }
              
            }
            
            # Format data for cumulative number of species
            
            species_acc_list_old <- vector("list")
            species_acc_list_sec <- vector("list")
            
            if(ncol(butterfly_summary_old) > 1){
              
              for (i in 1:ncol(butterfly_summary_old)){
                
                species_acc_list_old[[i]] <- butterfly_summary_old[,c(1:i)]
                species_acc_list_old[[1]] <- sum(species_acc_list_old[[1]]>0)
                
                species_acc_list_sec[[i]] <- butterfly_summary_sec[,c(1:i)]
                species_acc_list_sec[[1]] <- sum(species_acc_list_sec[[1]]>0)
                
              } 
              
              for (i in 2:length(species_acc_list_old)){
                
                species_acc_list_old[[i]] <- sum(rowSums(species_acc_list_old[[i]])>0)
                species_acc_list_sec[[i]] <- sum(rowSums(species_acc_list_sec[[i]])>0)
                
                
              }
              
            }
            
            # 5.12. Render species/diversity accumulation plots 
            
            if(ncol(butterfly_summary_old) > 1){
              
              # Species by samples
              
              species_by_sample_df_old <- data.frame(sample = seq(1, ncol(butterfly_summary_old), 1), 
                                                 count = unlist(species_acc_list_old))
              
              species_by_sample_df_old$forest_type <- "Old-growth"
              
              species_by_sample_df_sec <- data.frame(sample = seq(1, ncol(butterfly_summary_sec), 1), 
                                                 count = unlist(species_acc_list_sec))
              
              species_by_sample_df_sec$forest_type <- "Second-growth"
              
              species_by_sample_df <- rbind(species_by_sample_df_old, 
                                            species_by_sample_df_sec)
              
              output$species_by_sample_plot  <- 
                
                renderPlot({
                  
                  ggplot(species_by_sample_df,
                         aes(x = sample,
                             y = count, 
                             color = forest_type)) + 
                    geom_line(linewidth = 1) + 
                    geom_point(size = 4, 
                               shape = 21, 
                               fill = "white", 
                               stroke = 1) + 
                    theme_thomas()+
                    
                    xlab("Sample number")+
                    ylab("Species richness")+
                    
                    scale_x_continuous(expand = c(0.05, 0.5))+
                    scale_y_continuous(expand = c(0.1, 1))+
                    scale_color_manual(values = c("#769a6e", "#FFBF00"))
                  
                })
              
              
              # Species by individuals
              
              species_by_individuals_df_old <- data.frame(count = unlist(species_acc_list_old), 
                                                      individuals = sapply(individual_acc_list_old, function(x)
                                                        sum(unlist(x))))
              
              species_by_individuals_df_old$forest_type <- "Old-growth"
              
              species_by_individuals_df_sec <- data.frame(count = unlist(species_acc_list_sec), 
                                                          individuals = sapply(individual_acc_list_sec, function(x)
                                                            sum(unlist(x))))
              
              species_by_individuals_df_sec$forest_type <- "Second-growth"
              
              species_by_individuals_df <- rbind(species_by_individuals_df_old, 
                                                 species_by_individuals_df_sec)
              
              output$species_by_individuals_plot  <- 
                
                renderPlot({
                  
                  ggplot(species_by_individuals_df,
                         aes(x = individuals,
                             y = count, 
                             color = forest_type)) + 
                    geom_line(linewidth = 1) + 
                    geom_point(size = 4, 
                               shape = 21, 
                               fill = "white", 
                               stroke = 1) + 
                    theme_thomas()+
                    
                    xlab("Number of individuals")+
                    ylab("Species richness")+
                    
                    scale_x_continuous(expand = c(0.05, 0.5))+
                    scale_y_continuous(expand = c(0.1, 1))+
                    scale_color_manual(values = c("#769a6e", "#FFBF00"))
                  
                })
              
              # Shannon by samples
              
              shannon_by_samples_df_old <- data.frame(sample = seq(1, ncol(butterfly_summary_old), 1), 
                                                  count = sapply(individual_acc_list_old, function(x)
                                                    vegan::diversity(unlist(x), index = "shannon")))
              
              shannon_by_samples_df_old$forest_type <- "Old-growth"
              
              shannon_by_samples_df_sec <- data.frame(sample = seq(1, ncol(butterfly_summary_sec), 1), 
                                                      count = sapply(individual_acc_list_sec, function(x)
                                                        vegan::diversity(unlist(x), index = "shannon")))
              
              shannon_by_samples_df_sec$forest_type <- "Second-growth"
              
              shannon_by_samples_df <- rbind(shannon_by_samples_df_old, 
                                             shannon_by_samples_df_sec)
              
              output$shannon_by_samples_plot  <- 
                
                renderPlot({
                  
                  ggplot(shannon_by_samples_df,
                         aes(x = sample,
                             y = count, 
                             color = forest_type)) + 
                    geom_line(linewidth = 1) + 
                    geom_point(size = 4, 
                               shape = 21, 
                               fill = "white", 
                               stroke = 1) + 
                    theme_thomas()+
                    
                    xlab("Sample number")+
                    ylab("Shannon's H'")+
                    
                    scale_x_continuous(expand = c(0.05, 0.5))+
                    scale_y_continuous(expand = c(0.1, 1))+
                    scale_color_manual(values = c("#769a6e", "#FFBF00"))
                  
                })
              
              # Shannon by individuals
              
              shannon_by_individuals_df_old <- data.frame(individuals = sapply(individual_acc_list_old, function(x)
                sum(unlist(x))), 
                count = sapply(individual_acc_list_old, function(x)
                  vegan::diversity(unlist(x), index = "shannon")))
              
              shannon_by_individuals_df_old$forest_type <- "Old-growth"
              
              shannon_by_individuals_df_sec <- data.frame(individuals = sapply(individual_acc_list_sec, function(x)
                sum(unlist(x))), 
                count = sapply(individual_acc_list_sec, function(x)
                  vegan::diversity(unlist(x), index = "shannon")))
              
              shannon_by_individuals_df_sec$forest_type <- "Second-growth"
              
              shannon_by_individuals_df <- rbind(shannon_by_individuals_df_old, 
                                                 shannon_by_individuals_df_sec)
              
              output$shannon_by_individuals_plot  <- 
                
                renderPlot({
                  
                  ggplot(shannon_by_individuals_df,
                         aes(x = individuals,
                             y = count, 
                             color = forest_type)) + 
                    geom_line(linewidth = 1) + 
                    geom_point(size = 4, 
                               shape = 21, 
                               fill = "white", 
                               stroke = 1) + 
                    theme_thomas()+
                    
                    xlab("Number of individuals")+
                    ylab("Shannon's H'")+
                    
                    scale_x_continuous(expand = c(0.05, 0.5))+
                    scale_y_continuous(expand = c(0.1, 1))+
                    scale_color_manual(values = c("#769a6e", "#FFBF00"))
                  
                })
              
              # Simpson by samples
              
              simpson_by_samples_df_old <- data.frame(sample = seq(1, ncol(butterfly_summary_old), 1), 
                                                  count = sapply(individual_acc_list_old, function(x)
                                                    vegan::diversity(unlist(x), index = "simpson")))
              
              simpson_by_samples_df_old$forest_type <- "Old-growth"
              
              simpson_by_samples_df_sec <- data.frame(sample = seq(1, ncol(butterfly_summary_sec), 1), 
                                                  count = sapply(individual_acc_list_sec, function(x)
                                                    vegan::diversity(unlist(x), index = "simpson")))
              
              simpson_by_samples_df_sec$forest_type <- "Second-growth"
              
              simpson_by_samples_df <- rbind(simpson_by_samples_df_old, 
                                             simpson_by_samples_df_sec)
              
              output$simpson_by_samples_plot  <- 
                
                renderPlot({
                  
                  ggplot(simpson_by_samples_df,
                         aes(x = sample,
                             y = count, 
                             color = forest_type)) + 
                    geom_line(linewidth = 1) + 
                    geom_point(size = 4, 
                               shape = 21, 
                               fill = "white", 
                               stroke = 1) + 
                    theme_thomas()+
                    
                    xlab("Sample number")+
                    ylab("Simpson's D")+
                    
                    scale_x_continuous(expand = c(0.05, 0.5))+
                    scale_y_continuous(expand = c(0.1, 1))+
                    scale_color_manual(values = c("#769a6e", "#FFBF00"))
                  
                })
              
              # Simpson by individuals
              
              simpson_by_individuals_df_old <- data.frame(individuals = sapply(individual_acc_list_old, function(x)
                sum(unlist(x))), 
                count = sapply(individual_acc_list_old, function(x)
                  vegan::diversity(unlist(x), index = "simpson")))
              
              simpson_by_individuals_df_old$forest_type <- "Old-growth"
              
              simpson_by_individuals_df_sec <- data.frame(individuals = sapply(individual_acc_list_sec, function(x)
                sum(unlist(x))), 
                count = sapply(individual_acc_list_sec, function(x)
                  vegan::diversity(unlist(x), index = "simpson")))
              
              simpson_by_individuals_df_sec$forest_type <- "Second-growth"
              
              simpson_by_individuals_df <- rbind(simpson_by_individuals_df_old, 
                                                 simpson_by_individuals_df_sec)
              
              output$simpson_by_individuals_plot  <- 
                
                renderPlot({
                  
                  ggplot(simpson_by_individuals_df,
                         aes(x = individuals,
                             y = count, 
                             color = forest_type)) + 
                    geom_line(linewidth = 1) + 
                    geom_point(size = 4, 
                               shape = 21, 
                               fill = "white", 
                               stroke = 1) + 
                    theme_thomas()+
                    
                    xlab("Number of individuals")+
                    ylab("Simpson's D")+
                    
                    scale_x_continuous(expand = c(0.05, 0.5))+
                    scale_y_continuous(expand = c(0.1, 1))+
                    scale_color_manual(values = c("#769a6e", "#FFBF00"))
                  
                })
              
              
              # Evenness by samples
              
              evenness_by_samples_df_old <- data.frame(sample = seq(1, ncol(butterfly_summary_old), 1), 
                                                   count = sapply(individual_acc_list_old, function(x)
                                                     (vegan::diversity(x = unlist(x), index = "shannon")/log(sum(unlist(x>0))))))
              
              evenness_by_samples_df_old$forest_type <- "Old-growth"
              
              evenness_by_samples_df_sec <- data.frame(sample = seq(1, ncol(butterfly_summary_sec), 1), 
                                                       count = sapply(individual_acc_list_sec, function(x)
                                                         (vegan::diversity(x = unlist(x), index = "shannon")/log(sum(unlist(x>0))))))
              
              evenness_by_samples_df_sec$forest_type <- "Second-growth"
              
              evenness_by_samples_df <- rbind(evenness_by_samples_df_old, 
                                              evenness_by_samples_df_sec)
              
              
              output$evenness_by_samples_plot  <- 
                
                renderPlot({
                  
                  ggplot(evenness_by_samples_df,
                         aes(x = sample,
                             y = count, 
                             color = forest_type)) + 
                    geom_line(linewidth = 1) + 
                    geom_point(size = 4, 
                               shape = 21, 
                               fill = "white", 
                               stroke = 1) + 
                    theme_thomas()+
                    
                    xlab("Sample number")+
                    ylab("Pielou's Evenness")+
                    
                    scale_x_continuous(expand = c(0.05, 0.5))+
                    scale_y_continuous(expand = c(0.1, 1))+
                    scale_color_manual(values = c("#769a6e", "#FFBF00"))
                  
                })
              
              # Evenness by individuals
              
              evenness_by_individuals_df_old <- data.frame(individuals = sapply(individual_acc_list_old, function(x)
                sum(unlist(x))), 
                count = sapply(individual_acc_list_old, function(x)
                  (vegan::diversity(x = unlist(x), index = "shannon")/log(sum(unlist(x>0))))))
              
              evenness_by_individuals_df_old$forest_type <- "Old-growth"
              
              evenness_by_individuals_df_sec <- data.frame(individuals = sapply(individual_acc_list_sec, function(x)
                sum(unlist(x))), 
                count = sapply(individual_acc_list_sec, function(x)
                  (vegan::diversity(x = unlist(x), index = "shannon")/log(sum(unlist(x>0))))))
              
              evenness_by_individuals_df_sec$forest_type <- "Second-growth"
              
              evenness_by_individuals_df <- rbind(evenness_by_individuals_df_old, 
                                                  evenness_by_individuals_df_sec)
              
              
              output$evenness_by_individuals_plot <- 
                
                renderPlot({
                  
                  ggplot(evenness_by_individuals_df,
                         aes(x = individuals,
                             y = count, 
                             color = forest_type)) + 
                    geom_line(linewidth = 1) + 
                    geom_point(size = 4, 
                               shape = 21, 
                               fill = "white", 
                               stroke = 1) + 
                    theme_thomas()+
                    
                    xlab("Number of individuals")+
                    ylab("Pielou's Evenness")+
                    
                    scale_x_continuous(expand = c(0.05, 0.5))+
                    scale_y_continuous(expand = c(0.1, 1))+
                    scale_color_manual(values = c("#769a6e", "#FFBF00"))
                  
                })
              
              
              # Diversity (q = 0) by samples
              
              hill_q1_by_samples_df_old <- data.frame(sample = seq(1, ncol(butterfly_summary_old), 1),
                                                  count = sapply(individual_acc_list_old, function(x)
                                                    hillR::hill_taxa(comm = unlist(x), q = 1, MARGIN = 2)))
              
              hill_q1_by_samples_df_old$forest_type <- "Old-growth"
              
              hill_q1_by_samples_df_sec <- data.frame(sample = seq(1, ncol(butterfly_summary_sec), 1),
                                                      count = sapply(individual_acc_list_sec, function(x)
                                                        hillR::hill_taxa(comm = unlist(x), q = 1, MARGIN = 2)))
              
              hill_q1_by_samples_df_sec$forest_type <- "Second-growth"
              
              hill_q1_by_samples_df <- rbind(hill_q1_by_samples_df_old, 
                                             hill_q1_by_samples_df_sec)


              output$hill_q1_by_samples_plot  <-

                renderPlot({

                  ggplot(hill_q1_by_samples_df,
                         aes(x = sample,
                             y = count, 
                             color = forest_type)) +
                    geom_line(linewidth = 1) +
                    geom_point(size = 4,
                               shape = 21,
                               fill = "white",
                               stroke = 1) +
                    theme_thomas()+

                    xlab("Sample number")+
                    ylab("Diversity (q = 1)")+

                    scale_x_continuous(expand = c(0.05, 0.5))+
                    scale_y_continuous(expand = c(0.1, 1))+
                    scale_color_manual(values = c("#769a6e", "#FFBF00"))

                })
              
              # Diversity (q = 0) by individuals
              
              hill_q1_by_individuals_df_old <- data.frame(individuals = sapply(individual_acc_list_old, function(x)
                sum(unlist(x))),
                count = sapply(individual_acc_list_old, function(x)
                  hillR::hill_taxa(comm = unlist(x), q = 1, MARGIN = 2)))
              
              hill_q1_by_individuals_df_old$forest_type <- "Old-growth"
              
              hill_q1_by_individuals_df_sec <- data.frame(individuals = sapply(individual_acc_list_sec, function(x)
                sum(unlist(x))),
                count = sapply(individual_acc_list_sec, function(x)
                  hillR::hill_taxa(comm = unlist(x), q = 1, MARGIN = 2)))
              
              hill_q1_by_individuals_df_sec$forest_type <- "Second-growth"
              
              hill_q1_by_individuals_df <- rbind(hill_q1_by_individuals_df_old, 
                                                 hill_q1_by_individuals_df_sec)


              output$hill_q1_by_individuals_plot  <-

                renderPlot({

                  ggplot(hill_q1_by_individuals_df,
                         aes(x = individuals,
                             y = count, 
                             color = forest_type)) +
                    geom_line(linewidth = 1) +
                    geom_point(size = 4,
                               shape = 21,
                               fill = "white",
                               stroke = 1) +
                    theme_thomas()+

                    xlab("Number of individuals")+
                    ylab("Diversity (q = 1)")+

                    scale_x_continuous(expand = c(0.05, 0.5))+
                    scale_y_continuous(expand = c(0.1, 1))+
                    scale_color_manual(values = c("#769a6e", "#FFBF00"))

                })
              
              # Diversity (q = 2) by samples
              
              hill_q2_by_samples_df_old <- data.frame(sample = seq(1, ncol(butterfly_summary_old), 1),
                                                  count = sapply(individual_acc_list_old, function(x)
                                                    hillR::hill_taxa(comm = unlist(x), q = 2, MARGIN = 2)))
              
              hill_q2_by_samples_df_old$forest_type <- "Old-growth"
              
              hill_q2_by_samples_df_sec <- data.frame(sample = seq(1, ncol(butterfly_summary_sec), 1),
                                                      count = sapply(individual_acc_list_sec, function(x)
                                                        hillR::hill_taxa(comm = unlist(x), q = 2, MARGIN = 2)))
              
              hill_q2_by_samples_df_sec$forest_type <- "Second-growth"
              
              hill_q2_by_samples_df <- rbind(hill_q2_by_samples_df_old, 
                                             hill_q2_by_samples_df_sec)

              output$hill_q2_by_samples_plot  <-

                renderPlot({

                  ggplot(hill_q2_by_samples_df,
                         aes(x = sample,
                             y = count, 
                             color = forest_type)) +
                    geom_line(linewidth = 1) +
                    geom_point(size = 4,
                               shape = 21,
                               fill = "white",
                               stroke = 1) +
                    theme_thomas()+

                    xlab("Sample number")+
                    ylab("Diversity (q = 2)")+

                    scale_x_continuous(expand = c(0.05, 0.5))+
                    scale_y_continuous(expand = c(0.1, 1))+
                    scale_color_manual(values = c("#769a6e", "#FFBF00"))

                })
              
              # Diversity (q = 2) by individuals
              
              hill_q2_by_individuals_df_old <- data.frame(individuals = sapply(individual_acc_list_old, function(x)
                sum(unlist(x))),
                count = sapply(individual_acc_list_old, function(x)
                  hillR::hill_taxa(comm = unlist(x), q = 2, MARGIN = 2)))
              
              hill_q2_by_individuals_df_old$forest_type <- "Old-growth"
              
              hill_q2_by_individuals_df_sec <- data.frame(individuals = sapply(individual_acc_list_sec, function(x)
                sum(unlist(x))),
                count = sapply(individual_acc_list_sec, function(x)
                  hillR::hill_taxa(comm = unlist(x), q = 2, MARGIN = 2)))
              
              hill_q2_by_individuals_df_sec$forest_type <- "Second-growth"
              
              hill_q2_by_individuals_df <- rbind(hill_q2_by_individuals_df_old, 
                                                 hill_q2_by_individuals_df_sec)

              output$hill_q2_by_individuals_plot  <-

                renderPlot({

                  ggplot(hill_q2_by_individuals_df,
                         aes(x = individuals,
                             y = count, 
                             color = forest_type)) +
                    geom_line(linewidth = 1) +
                    geom_point(size = 4,
                               shape = 21,
                               fill = "white",
                               stroke = 1) +
                    theme_thomas()+

                    xlab("Number of individuals")+
                    ylab("Diversity (q = 2)")+

                    scale_x_continuous(expand = c(0.05, 0.5))+
                    scale_y_continuous(expand = c(0.1, 1))+
                    scale_color_manual(values = c("#769a6e", "#FFBF00"))

                })
              
              # Evenness (Hill-based) by samples
              
              hill_evenness_by_samples_df_old <- data.frame(sample = seq(1, ncol(butterfly_summary_old), 1),
                                                        count = sapply(individual_acc_list_old, function(x)
                                                          ((hillR::hill_taxa(comm = unlist(x), q = 2, MARGIN = 2)-1)/
                                                             (hillR::hill_taxa(comm = unlist(x), q = 0, MARGIN = 2)-1))))
              
              hill_evenness_by_samples_df_old$forest_type <- "Old-growth"
              
              hill_evenness_by_samples_df_sec <- data.frame(sample = seq(1, ncol(butterfly_summary_sec), 1),
                                                            count = sapply(individual_acc_list_sec, function(x)
                                                              ((hillR::hill_taxa(comm = unlist(x), q = 2, MARGIN = 2)-1)/
                                                                 (hillR::hill_taxa(comm = unlist(x), q = 0, MARGIN = 2)-1))))
              
              hill_evenness_by_samples_df_sec$forest_type <- "Second-growth"
              
              hill_evenness_by_samples_df <- rbind(hill_evenness_by_samples_df_old, 
                                                   hill_evenness_by_samples_df_sec)

              output$hill_evenness_by_samples_plot  <-

                renderPlot({

                  ggplot(hill_evenness_by_samples_df,
                         aes(x = sample,
                             y = count, 
                             color = forest_type)) +
                    geom_line(linewidth = 1) +
                    geom_point(size = 4,
                               shape = 21,
                               fill = "white",
                               stroke = 1) +
                    theme_thomas()+

                    xlab("Sample number")+
                    ylab("Evenness (Hill-based)")+

                    scale_x_continuous(expand = c(0.05, 0.5))+
                    scale_y_continuous(expand = c(0.1, 1))+
                    scale_color_manual(values = c("#769a6e", "#FFBF00"))

                })
              
              # Evenness (Hill-based) by individuals
              
              hill_evenness_by_individuals_df_old <- data.frame(individuals = sapply(individual_acc_list_old, function(x)
                sum(unlist(x))),
                count = sapply(individual_acc_list_old, function(x)
                  ((hillR::hill_taxa(comm = unlist(x), q = 2, MARGIN = 2)-1)/
                     (hillR::hill_taxa(comm = unlist(x), q = 0, MARGIN = 2)-1))))
              
              hill_evenness_by_individuals_df_old$forest_type <- "Old-growth"
              
              hill_evenness_by_individuals_df_sec <- data.frame(individuals = sapply(individual_acc_list_sec, function(x)
                sum(unlist(x))),
                count = sapply(individual_acc_list_sec, function(x)
                  ((hillR::hill_taxa(comm = unlist(x), q = 2, MARGIN = 2)-1)/
                     (hillR::hill_taxa(comm = unlist(x), q = 0, MARGIN = 2)-1))))
              
              hill_evenness_by_individuals_df_sec$forest_type <- "Second-growth"
              
              hill_evenness_by_individuals_df <- rbind(hill_evenness_by_individuals_df_old, 
                                                       hill_evenness_by_individuals_df_sec)

              output$hill_evenness_by_individuals_plot  <-

                renderPlot({

                  ggplot(hill_evenness_by_individuals_df,
                         aes(x = individuals,
                             y = count, 
                             color = forest_type)) +
                    geom_line(linewidth = 1) +
                    geom_point(size = 4,
                               shape = 21,
                               fill = "white",
                               stroke = 1) +
                    theme_thomas()+

                    xlab("Number of individuals")+
                    ylab("Evenness (Hill-based)")+

                    scale_x_continuous(expand = c(0.05, 0.5))+
                    scale_y_continuous(expand = c(0.1, 1))+
                    scale_color_manual(values = c("#769a6e", "#FFBF00"))

                })
              
              # 5.13. Calculate the traditional diversity metrics
                
                sample_abundance_list_old <- vector("list")
                
                for (i in 1:ncol(butterfly_dataset_old())){
                  
                  sample_abundance_list_old[[i]] <- butterfly_dataset_old()[,c(1:i)]
                  
                } 
                
                for (i in 2:length(species_acc_list_old)){
                  
                  sample_abundance_list_old[[i]] <- rowSums(sample_abundance_list_old[[i]])
                  
                }
                
                sample_abundance_list_sec <- vector("list")
                
                for (i in 1:ncol(butterfly_dataset_sec())){
                  
                  sample_abundance_list_sec[[i]] <- butterfly_dataset_sec()[,c(1:i)]
                  
                } 
                
                for (i in 2:length(species_acc_list_sec)){
                  
                  sample_abundance_list_sec[[i]] <- rowSums(sample_abundance_list_sec[[i]])
                  
                }
                
                # Species richness
                
                species_richness_old <- lapply(sample_abundance_list_old, 
                                           function(x) sum(x>0))
                
                species_richness_sec <- lapply(sample_abundance_list_sec, 
                                               function(x) sum(x>0))
                
                # Shannon index
                
                shannon_index_old <- lapply(sample_abundance_list_old, 
                                        function(x) round(vegan::diversity(x = x, index = "shannon"),digits = 3))
                
                shannon_index_sec <- lapply(sample_abundance_list_sec, 
                                            function(x) round(vegan::diversity(x = x, index = "shannon"),digits = 3))
                
                
                # Simpson index
                
                simpson_index_old <- lapply(sample_abundance_list_old, 
                                        function(x) round(vegan::diversity(x = x, index = "simpson"), digits = 3))
                
                simpson_index_sec <- lapply(sample_abundance_list_sec, 
                                            function(x) round(vegan::diversity(x = x, index = "simpson"), digits = 3))
                
                
                # Species evenness
                
                species_evenness_old <- lapply(sample_abundance_list_old, 
                                           function(x) round((vegan::diversity(x = x, index = "shannon")/log(sum(x>0))), digits = 3))
                
                species_evenness_sec <- lapply(sample_abundance_list_sec, 
                                               function(x) round((vegan::diversity(x = x, index = "shannon")/log(sum(x>0))), digits = 3))
                
                # q = 0
                
                hill_q0_old <- lapply(sample_abundance_list_old, 
                                  function(x) round(hillR::hill_taxa(comm = x, 
                                                                     q = 0, 
                                                                     MARGIN = 2), digits = 3))
                
                hill_q0_sec <- lapply(sample_abundance_list_sec, 
                                  function(x) round(hillR::hill_taxa(comm = x, 
                                                                     q = 0, 
                                                                     MARGIN = 2), digits = 3))
                
                # q = 1
                
                hill_q1_old <- lapply(sample_abundance_list_old, 
                                  function(x) round(hillR::hill_taxa(comm = x, 
                                                                     q = 1, 
                                                                     MARGIN = 2), digits = 3))
                
                hill_q1_sec <- lapply(sample_abundance_list_sec, 
                                      function(x) round(hillR::hill_taxa(comm = x, 
                                                                         q = 1, 
                                                                         MARGIN = 2), digits = 3))
                
                
                # q = 2
                
                hill_q2_old <- lapply(sample_abundance_list_old, 
                                  function(x) round(hillR::hill_taxa(comm = x,
                                                                     q = 2, 
                                                                     MARGIN = 2), digits = 3))
                
                hill_q2_sec <- lapply(sample_abundance_list_sec, 
                                      function(x) round(hillR::hill_taxa(comm = x,
                                                                         q = 2, 
                                                                         MARGIN = 2), digits = 3))
                
                # Hill-based evenness
                
                hill_evenness_old <- lapply(sample_abundance_list_old, 
                                        function(x) round(((hillR::hill_taxa(comm = x, 
                                                                             q = 2, 
                                                                             MARGIN = 2)-1)) / (hillR::hill_taxa(comm = x, 
                                                                                                                 q = 0, 
                                                                                                                 MARGIN = 2)-1), 
                                                          digits = 3))
                
                hill_evenness_sec <- lapply(sample_abundance_list_sec, 
                                            function(x) round(((hillR::hill_taxa(comm = x, 
                                                                                 q = 2, 
                                                                                 MARGIN = 2)-1)) / (hillR::hill_taxa(comm = x, 
                                                                                                                     q = 0, 
                                                                                                                     MARGIN = 2)-1), 
                                                              digits = 3))
                
                
                traditional_table_old <- rbind(unlist(species_richness_old), 
                                           unlist(shannon_index_old), 
                                           unlist(simpson_index_old), 
                                           unlist(species_evenness_old))
                
                traditional_table_sec <- rbind(unlist(species_richness_sec), 
                                               unlist(shannon_index_sec), 
                                               unlist(simpson_index_sec), 
                                               unlist(species_evenness_sec))
                
                rownames(traditional_table_old) <- c("Species richness", 
                                                 "Shannon's H'", 
                                                 "Simpson's D", 
                                                 "Pielou's evenness")
                
                rownames(traditional_table_sec) <- c("Species richness", 
                                                     "Shannon's H'", 
                                                     "Simpson's D", 
                                                     "Pielou's evenness")
                
                colnames(traditional_table_old) <- paste0("Sample_", 
                                                      seq(1, ncol(butterfly_dataset_old()), 1))
                
                colnames(traditional_table_sec) <- paste0("Sample_", 
                                                          seq(1, ncol(butterfly_dataset_sec()), 1))
                
                hill_table_old <- rbind(unlist(hill_q0_old), 
                                    unlist(hill_q1_old), 
                                    unlist(hill_q2_old), 
                                    unlist(hill_evenness_old))
                
                hill_table_sec <- rbind(unlist(hill_q0_sec), 
                                        unlist(hill_q1_sec), 
                                        unlist(hill_q2_sec), 
                                        unlist(hill_evenness_sec))
                
                rownames(hill_table_old) <- c("Diversity (q = 0)", 
                                          "Diversity (q = 1)", 
                                          "Diversity (q = 2)", 
                                          "Hill evenness")
                
                rownames(hill_table_sec) <- c("Diversity (q = 0)", 
                                              "Diversity (q = 1)", 
                                              "Diversity (q = 2)", 
                                              "Hill evenness")
                
                colnames(hill_table_old) <- paste0("Sample_", 
                                               seq(1, ncol(butterfly_dataset_old()), 1))
                
                colnames(hill_table_sec) <- paste0("Sample_", 
                                                   seq(1, ncol(butterfly_dataset_sec()), 1))
                
                
                output$traditional_table <- renderDT({
                  datatable(
                    traditional_table_old,
                    rownames = TRUE,
                    options = list(scrollX = TRUE, scrollY = "250px")
                  )
                })
                
                output$traditional_table_2 <- renderDT({
                  datatable(
                    traditional_table_sec,
                    rownames = TRUE,
                    options = list(scrollX = TRUE, scrollY = "250px")
                  )
                })
                
                
                output$hill_table <- renderDT({
                  datatable(
                    hill_table_old,
                    rownames = TRUE,
                    options = list(scrollX = TRUE, scrollY = "250px")
                  )
                })
                
                output$hill_table_2 <- renderDT({
                  datatable(
                    hill_table_sec,
                    rownames = TRUE,
                    options = list(scrollX = TRUE, scrollY = "250px")
                  )
                })
              
            }
            
            
          }
          
        }
        
    })
    
    # 5.13. Prepare individual-based rarefaction data
    
    observeEvent(input$rarefaction_button, {
      
      if(input$forest_type == "Old-Growth" | input$forest_type == "Second-Growth"){
      
      # Individual-based rarefaction plot
      
      rar_abun_data <- vector("list")
      rar_abun_data[[1]] <- as.vector(rowSums(butterfly_dataset()))
      
      names(rar_abun_data) <- as.character(input$forest_type)
      
      iNEXT_individuals <- iNEXT::iNEXT(x = rar_abun_data, 
                                        q = c(0, 1, 2), 
                                        datatype = "abundance")
      
      iNEXT_individuals_siz_table <- iNEXT_individuals[[2]]$size_based
      iNEXT_individuals_cov_table <- iNEXT_individuals[[2]]$coverage_based

      output$rar_table_abun_siz <-

        renderDT({

          datatable(
            iNEXT_individuals_siz_table,
            rownames = TRUE,
            options = list(scrollX = TRUE, scrollY = "250px"))

        })

      output$rar_table_abun_cov <-

        renderDT({

          datatable(
            iNEXT_individuals_cov_table,
            rownames = TRUE,
            options = list(scrollX = TRUE, scrollY = "250px"))

        })
      
      output$rar_plot_abun_siz  <- 
        
        renderPlot({
          
          iNEXT::ggiNEXT(iNEXT_individuals, 
                         type=1, 
                         facet.var = "Order.q") +
            theme_thomas() + 
            
            theme(legend.position = "none")
          
          
        })
      
      output$rar_plot_abun_cov <-
        
        renderPlot({
          
          iNEXT::ggiNEXT(iNEXT_individuals,
                         type=3,
                         facet.var="Order.q") +
            
            theme_thomas() +
            
            theme(legend.position = "none")
          
          
        })
      
      # Sample-based rarefaction plot
      
      incidence_data <- butterfly_dataset()
      incidence_data[incidence_data > 1] <- 1
      
      nsamples <- ncol(incidence_data)
      ndetections <- rowSums(incidence_data)
      incidence_data <- c(nsamples, ndetections)
      
      incidence_data <- list(forest = incidence_data)
      names(incidence_data) <- as.character(input$forest_type)
      
      iNEXT_samples <- iNEXT::iNEXT(x = incidence_data, 
                                    q = c(0, 1, 2), 
                                    datatype = "incidence_freq")
      
      iNEXT_samples_siz_table <- iNEXT_samples[[2]]$size_based
      iNEXT_samples_cov_table <- iNEXT_samples[[2]]$coverage_based

      output$rar_table_inc_siz <-

        renderDT({

          datatable(
            iNEXT_samples_siz_table,
            rownames = TRUE,
            options = list(scrollX = TRUE, scrollY = "250px"))

        })

      output$rar_table_inc_cov <-

        renderDT({

          datatable(
            iNEXT_samples_cov_table,
            rownames = TRUE,
            options = list(scrollX = TRUE, scrollY = "250px"))

        })
      
      output$rar_plot_inc_siz  <- 
        
        renderPlot({
          
          iNEXT::ggiNEXT(iNEXT_samples, 
                         type=1, 
                         facet.var = "Order.q") +
            theme_thomas() + 
            
            theme(legend.position = "none")
          
          
        })
      
      output$rar_plot_inc_cov <-
        
        renderPlot({
          
          iNEXT::ggiNEXT(iNEXT_samples,
                         type=3,
                         facet.var="Order.q") +
            
            theme_thomas() +
            
            theme(legend.position = "none")
          
          
        })
      
      } else{
        
        if(input$forest_type == "Both"){
          
          # 5.13. Prepare rarefaction plots
          
          
          
          # Individual-based rarefaction plot
          
          rar_abun_data_old <- vector("list")
          rar_abun_data_old[[1]] <- as.vector(rowSums(butterfly_dataset_old()))
          names(rar_abun_data_old) <- "Old-growth"
          
          rar_abun_data_sec <- vector("list")
          rar_abun_data_sec[[1]] <- as.vector(rowSums(butterfly_dataset_sec()))
          names(rar_abun_data_sec) <- "Second-growth"
          
          rar_abun_data <- vector("list")
          rar_abun_data[[1]] <- unlist(rar_abun_data_old[[1]])
          rar_abun_data[[2]] <- unlist(rar_abun_data_sec[[1]])
          names(rar_abun_data) <- c("Old-growth", "Second-growth")
          
          iNEXT_individuals <- iNEXT::iNEXT(x = rar_abun_data, 
                                            q = c(0, 1, 2), 
                                            datatype = "abundance")
          
          iNEXT_individuals_siz_table <- iNEXT_individuals[[2]]$size_based
          iNEXT_individuals_cov_table <- iNEXT_individuals[[2]]$coverage_based
          
          output$rar_table_abun_siz <-
            
            renderDT({
              
              datatable(
                iNEXT_individuals_siz_table,
                rownames = TRUE,
                options = list(scrollX = TRUE, scrollY = "250px"))
              
            })
          
          output$rar_table_abun_cov <-
            
            renderDT({
              
              datatable(
                iNEXT_individuals_cov_table,
                rownames = TRUE,
                options = list(scrollX = TRUE, scrollY = "250px"))
              
            })
          
          output$rar_plot_abun_siz  <- 
            
            renderPlot({
              
              iNEXT::ggiNEXT(iNEXT_individuals, 
                             type=1, 
                             facet.var = "Order.q", 
                             color.var = "Assemblage") +
                theme_thomas() + 
                
                theme(legend.position = "none")+
                scale_color_manual(values = c("#769a6e", "#FFBF00"))+
                scale_fill_manual(values = c("#769a6e", "#FFBF00"))
              
              
            })
          
          output$rar_plot_abun_cov <-
            
            renderPlot({
              
              iNEXT::ggiNEXT(iNEXT_individuals,
                             type=3,
                             facet.var="Order.q",
                             color.var = "Assemblage") +
                
                theme_thomas() +
                
                theme(legend.position = "none")+
                scale_color_manual(values = c("#769a6e", "#FFBF00"))+
                scale_fill_manual(values = c("#769a6e", "#FFBF00"))
              
              
            })
          
          # Sample-based rarefaction plot
          
          incidence_data_old <- butterfly_dataset_old()
          incidence_data_old[incidence_data_old > 1] <- 1
          
          nsamples_old <- ncol(incidence_data_old)
          ndetections_old <- rowSums(incidence_data_old)
          incidence_data_old <- c(nsamples_old, ndetections_old)
          
          incidence_data_sec <- butterfly_dataset_sec()
          incidence_data_sec[incidence_data_sec > 1] <- 1
          
          nsamples_sec <- ncol(incidence_data_sec)
          ndetections_sec <- rowSums(incidence_data_sec)
          incidence_data_sec <- c(nsamples_sec, ndetections_sec)
          
          incidence_data <- vector("list")
          incidence_data[[1]] <- unlist(incidence_data_old)
          incidence_data[[2]] <- unlist(incidence_data_sec)
          names(incidence_data) <- c("Old-growth", "Second-growth")
          
          
          iNEXT_samples <- iNEXT::iNEXT(x = incidence_data,
                                        q = c(0, 1, 2),
                                        datatype = "incidence_freq")
          
          iNEXT_samples_siz_table <- iNEXT_samples[[2]]$size_based
          iNEXT_samples_cov_table <- iNEXT_samples[[2]]$coverage_based
          
          output$rar_table_inc_siz <-
            
            renderDT({
              
              datatable(
                iNEXT_samples_siz_table,
                rownames = TRUE,
                options = list(scrollX = TRUE, scrollY = "250px"))
              
            })
          
          output$rar_table_inc_cov <-
            
            renderDT({
              
              datatable(
                iNEXT_samples_cov_table,
                rownames = TRUE,
                options = list(scrollX = TRUE, scrollY = "250px"))
              
            })
          
          output$rar_plot_inc_siz  <-
            
            renderPlot({
              
              iNEXT::ggiNEXT(iNEXT_samples,
                             type=1,
                             facet.var = "Order.q", 
                             color.var = "Assemblage") +
                theme_thomas() +
                
                theme(legend.position = "none")+
                scale_color_manual(values = c("#769a6e", "#FFBF00"))+
                scale_fill_manual(values = c("#769a6e", "#FFBF00"))
              
              
            })
          
          output$rar_plot_inc_cov <-
            
            renderPlot({
              
              iNEXT::ggiNEXT(iNEXT_samples,
                             type=3,
                             facet.var="Order.q", 
                             color.var = "Assemblage") +
                
                theme_thomas() +
                
                theme(legend.position = "none")+
                scale_color_manual(values = c("#769a6e", "#FFBF00"))+
                scale_fill_manual(values = c("#769a6e", "#FFBF00"))
              
              
            })
          
          
        }
        
        
      }
      
    })
}

# Run the app
shinyApp(ui, server)