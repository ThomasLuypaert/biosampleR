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
num_individuals_old_growth <- 10000
num_individuals_second_growth <- 6500

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
abundance_old_growth[abundance_old_growth == 0] <- round(runif(sum(abundance_old_growth == 0), min = 1, max = 90))

# Assign abundance for common species in second-growth
common_species_second_growth <- sample(second_growth_community, num_common_species_second_growth, replace = FALSE)
abundance_second_growth[second_growth_community %in% common_species_second_growth] <- round(runif(num_common_species_second_growth, min = 100, max = 500))

abundance_second_growth[abundance_second_growth == 0] <- round(runif(sum(abundance_second_growth == 0), min = 1, max = 50))

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
                         selectInput("forest_type", "Select Forest Type", c("Old-Growth", "Second-Growth")),
                         actionButton("sample_button", "Draw Random Sample"),
                         actionButton("reset_button", "Reset"),
                         br(), 
                         br(), 
                         textOutput("counter_text"), 
                         textOutput("sampled_text")
                         
                     ),
                     mainPanel(plotOutput("butterfly_plot_1"),
                               DTOutput("sample_table")
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
                    br(),
                    DTOutput("traditional_table"), 
                    br(), 
                    br(),
                    h4(strong("B. Hill-based diversity metrics")),
                    br(), 
                    br(),
                    DTOutput("hill_table"),
                    br(), 
                ),
                tabPanel("Sampling procedure",
                    h3("Explore sampling procedure"),
                    br(),
                    
                    h4(strong("A. Species-by-sample")),
                    br(),
                    plotOutput("species_sample_plot"), 
                    br(),
                    plotOutput("species_acc_plot"), 
                    br(), 
                    h4(strong("B. Diversity-by-sample")),
                    br()
                    # ... (add UI elements specific to sampling effort)
                ),
                tabPanel("Abundance Distribution",
                    h3("Explore Abundance Distribution"),
                    plotOutput("butterfly_plot_2")
                    # ... (add UI elements specific to abundance distribution)
                ),
                tabPanel("Rarefaction",
                         h3("Explore rarefaction results"),
                         plotOutput("rar_plot_abun_siz")
                ),
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
        
        # Update the UI element displaying the sampled_text
        output$sampled_text <- renderText({
            paste("Number of individuals sampled:", 0)
        })
        
        # Clear the plot, table, and butterfly_data
        output$butterfly_plot_1 <- renderPlot(NULL)
        output$butterfly_plot_2 <- renderPlot(NULL)
        output$sample_table <- renderDT(NULL)
        butterfly_data$data <- vector("list")
        output$traditional_table <- renderDT(NULL)
        output$hill_table <- renderDT(NULL)
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
                       "Second-Growth" = dataset_second_growth)
            )
            
            # Update the UI element displaying the sampled_text
            output$sampled_text <- renderText({
                paste("Number of individuals sampled:", 0)
            })
            
            # Clear the plot and table
            output$butterfly_plot_1 <- renderPlot(NULL)
            output$butterfly_plot_2 <- renderPlot(NULL)
            output$sample_table <- renderDT(NULL)
            butterfly_data$data <- vector("list")
            
            output$traditional_table <- renderDT(NULL)
            output$hill_table <- renderDT(NULL)
        }
    }, priority = 1000)
    
    # 3. Create reactive object to save sampling data
    
    butterfly_data <- reactiveValues(data = list())
    
    # 4. Create reactive object to store community data after sampling
    
    selected_dataset <- reactiveVal(NULL)
    
    observeEvent(input$forest_type, {
        selected_dataset(
            switch(input$forest_type,
                   "Old-Growth" = dataset_old_growth,
                   "Second-Growth" = dataset_second_growth)
        )
    })
    
    # 5. Create events when sample button is clicked 
    
    
    observeEvent(input$sample_button, {
        
        # 5.1. Determine how many individuals will be sampled
        
        num_individuals_sample <- sample(c(1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 8, 8, 9, 9, 10, 11, 12), 
                                         size = 1)
        
        
        # 5.2. Sample the community
        
        sample_result <- sample(selected_dataset(), size = num_individuals_sample, replace = FALSE)
        
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
        
        # 5.6. Create counter of how many individuals sampled
        
        output$sampled_text <- renderText({
            paste("Number of individuals:", sum(unlist(butterfly_summary)))
        })
        
        
        # 5.6. Display sample data table
        output$sample_table <- renderDT({
            datatable(
                butterfly_summary,
                rownames = TRUE,
                options = list(scrollX = TRUE, scrollY = "400px")
            )
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
                geom_line(color = "#769a6e",  
                          linewidth = 1) + 
                geom_point(size = 4, 
                           shape = 21, 
                           color = "#769a6e", 
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
        
        sample_plot_data <- data.frame(count = colSums(butterfly_summary), 
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
              ggtitle("Number of species per sample by sampling effort")+
              
              scale_x_continuous(expand = c(0.05, 0.5))+
              scale_y_continuous(expand = c(0.1, 1))
            
          })
        
        
        # 5.11. Format data for species accumulation plot
        
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
        
      
        # 5.12. Render species-by-sample plot
        
        if(ncol(butterfly_summary) > 1){
          
         species_acc_df <- data.frame(sample = seq(1, ncol(butterfly_summary), 1), 
                                       count = unlist(species_acc_list))
        
        output$species_acc_plot  <- 
          
          renderPlot({
            
            ggplot(species_acc_df,
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
              ylab("Total species richness")+
              ggtitle("Cumulative species richness by sampling effort")+
              
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
        
        
        # 5.13. Prepare individual-based rarefaction data
        
        if(ncol(butterfly_summary) > 3){
          
          rar_abun_data <- vector("list")
          rar_abun_data[[1]] <- as.vector(rowSums(butterfly_summary))
          
          names(rar_abun_data) <- "data"
          
          iNEXT <- iNEXT::iNEXT(x = rar_abun_data, 
                                q = c(0, 1, 2), 
                                datatype = "abundance")
          
          output$rar_plot_abun_siz  <- 
            
            renderPlot({
              
              iNEXT::ggiNEXT(iNEXT, 
                             type=1, 
                             facet.var = "Order.q") +
                theme_thomas() + 
                
                theme(legend.position = "none")
            
              
            })
        
        }
        
    })
}

# Run the app
shinyApp(ui, server)