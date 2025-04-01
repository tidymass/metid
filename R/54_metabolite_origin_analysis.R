# ####example
# library(r4projects)
# setwd(r4projects::get_project_wd())
# setwd("other_files/metabolite_origin_analysis/")
# # rm(list = ls())
# # gc()
# load("pregnancy_urine_metabolomics.rda")
# load("hmdb_ms2.rda")
#
# # load("ms1_database.rda")
# # annotation_table <-
# #   pregnancy_urine_metabolomics@variable_info
# #
# # annotation_table <-
# #   annotation_table %>%
# #   dplyr::filter(!is.na(HMDB.ID)) %>%
# #   dplyr::distinct(HMDB.ID, .keep_all = TRUE) %>%
# #   dplyr::left_join(
# #     ms1_database %>%
# #       dplyr::filter(!is.na(HMDB_ID)) %>%
# #       dplyr::distinct(HMDB_ID, .keep_all = TRUE) %>%
# #       dplyr::select(c(HMDB_ID, from_human:from_which_food)),
# #     by = c("HMDB.ID" = "HMDB_ID")
# #   )
# #
# # save(annotation_table, file = "annotation_table.rda")
# load("annotation_table.rda")
#
# pregnancy_urine_metabolomics@annotation_table <-
#   annotation_table
#
# object <-
#   pregnancy_urine_metabolomics
#
# object <-
#   analyze_metabolite_origins(object = object)
#
# metabolite_origin_upsetplot(object = object)
#
#
# ###metabolite_origin_network
# metabolite_origin_network(
#   object = object,
#   metabolite_id = object@annotation_table$Lab.ID[c(421)],
#   top_specific_source = 3
# )
#
# metabolite_origin_network(
#   object = object,
#   metabolite_id = object@annotation_table$Lab.ID[c(421, 422)],
#   top_specific_source = 3
# )
#
# source_network(
#   object = object,
#   source_id = c("Food"),
#   top_specific_source = 3
# )
#
# specific_source_network(
#   object = object,
#   specific_source_id = c("Urine"),
#   top_specific_source = 3
# )


#' Filter a dataset to include only metabolites of known origin
#'
#' @description
#' This function filters a mass dataset to keep only metabolites with known
#' origin information, removing entries where the 'from_human' field is NA.
#'
#' @param object A data.frame or mass_dataset object containing metabolite annotation data
#'
#' @return The filtered object with the same class as the input
#'
#' @author Xiaotao Shen \email{xiaotao.shen@@outlook.com}
#'
#' @importFrom dplyr distinct filter
#'
#' @export

analyze_metabolite_origins <-
  function(object) {
    check_object4metablite_origin(object)
    
    if (is(object, "mass_dataset")) {
      annotation_table <-
        object@annotation_table
    }
    
    annotation_table <-
      annotation_table %>%
      dplyr::distinct(Lab.ID, .keep_all = TRUE) %>%
      dplyr::filter(!is.na(from_human))
    
    object@annotation_table <-
      annotation_table
    return(object)
  }


#' Create a network visualization of metabolite sources
#'
#' @description
#' Generates a network plot visualizing the relationships between metabolites and their sources.
#' The network shows metabolites connected to their origins (human, bacteria, plant, etc.).
#'
#' @param object A mass_dataset object containing metabolite annotation data
#' @param circle Logical specifying whether to display the network in a circular layout (default: TRUE)
#'
#' @return A ggplot2 object representing the source network visualization
#'
#' @author Xiaotao Shen \email{xiaotao.shen@@outlook.com}
#'
#' @importFrom ggraph ggraph geom_edge_diagonal geom_node_point geom_node_text scale_edge_color_manual theme_graph create_layout
#' @importFrom tidygraph tbl_graph activate mutate
#' @importFrom igraph bipartite_mapping V
#' @importFrom dplyr select filter mutate arrange distinct rename pull left_join
#' @importFrom stringr str_split
#' @importFrom purrr map
#' @importFrom ggplot2 scale_fill_manual theme element_blank
#' @importFrom scales scale_size_continuous
#'
#' @export

source_metabolite_network <-
  function(object, circle = TRUE) {
    check_object4metablite_origin(object)
    temp_data <-
      object@annotation_table %>%
      dplyr::select(
        Lab.ID,
        from_human,
        from_which_part,
        from_bacteria,
        from_which_bacteria,
        from_plant,
        from_which_plant,
        from_animal,
        from_which_animal,
        from_environment,
        from_which_environment,
        from_drug,
        from_which_drug,
        from_food,
        from_which_food
      )
    
    colnames(temp_data) <-
      c(
        "Lab.ID",
        "Human",
        "Human_name",
        "Bacteria",
        "Bacteria_name",
        "Plant",
        "Plant_name",
        "Animal",
        "Animal_name",
        "Environment",
        "Environment_name",
        "Drug",
        "Drug_name",
        "Food",
        "Food_name"
      )
    
    temp_data$Human_name[temp_data$Human_name == "Unknown"] <-
      paste("Human", "Unknown", sep = "_")
    
    temp_data$Bacteria_name[temp_data$Bacteria_name == "Unknown"] <-
      paste("Bacteria", "Unknown", sep = "_")
    
    temp_data$Plant_name[temp_data$Plant_name == "Unknown"] <-
      paste("Plant", "Unknown", sep = "_")
    
    temp_data$Animal_name[temp_data$Animal_name == "Unknown"] <-
      paste("Animal", "Unknown", sep = "_")
    
    temp_data$Environment_name[temp_data$Environment_name == "Unknown"] <-
      paste("Environment", "Unknown", sep = "_")
    
    temp_data$Drug_name[temp_data$Drug_name == "Unknown"] <-
      paste("Drug", "Unknown", sep = "_")
    
    temp_data$Food_name[temp_data$Food_name == "Unknown"] <-
      paste("Food", "Unknown", sep = "_")
    
    Source <-
      unique(temp_data$Lab.ID) %>%
      purrr::map(function(x) {
        temp <-
          temp_data %>%
          dplyr::filter(Lab.ID == x)
        colnames(temp)[which(as.character(temp) == "Yes")]
      })
    
    temp_data2 <-
      temp_data %>%
      dplyr::select(Lab.ID)
    
    temp_data2$source <- Source
    
    ###network
    ###source information
    temp_data3 <-
      seq_len(nrow(temp_data2)) %>%
      purrr::map(function(i) {
        temp <-
          tryCatch(
            data.frame(Lab.ID = temp_data2$Lab.ID[i], Source = temp_data2$source[[i]]),
            error = function(e) {
              NULL
            }
          )
        temp
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    #####specific source
    temp_data4 <-
      unique(temp_data3$Lab.ID) %>%
      purrr::map(function(x) {
        temp <-
          temp_data3 %>%
          dplyr::filter(Lab.ID == x)
        index <-
          temp %>%
          dplyr::pull(Source) %>%
          paste0("_name")
        
        specific_source <-
          temp_data %>%
          dplyr::filter(Lab.ID == x) %>%
          dplyr::select(all_of(index)) %>%
          as.character() %>%
          purrr::map(function(y) {
            temp <-
              stringr::str_split(y, "\\{\\}")[[1]]
            temp
          })
        
        temp <-
          seq_len(nrow(temp)) %>%
          purrr::map(function(i) {
            data.frame(
              Lab.ID = temp$Lab.ID[i],
              Source = temp$Source[i],
              Specific_source = specific_source[[i]]
            )
          }) %>%
          do.call(rbind, .) %>%
          as.data.frame()
        
        temp
        
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    edge_data <-
      temp_data3 %>%
      dplyr::rename(from = Lab.ID, to = Source) %>%
      dplyr::mutate(edge_class = to)
    
    node_data <-
      data.frame(name = unique(edge_data$from)) %>%
      dplyr::left_join(object@annotation_table[, c("Lab.ID", "Compound.name")], by = c("name" = "Lab.ID")) %>%
      dplyr::mutate(class = "metabolite")
    
    node_data <-
      rbind(node_data, data.frame(
        name = c(unique(edge_data$to)),
        Compound.name = c(unique(edge_data$to)),
        class = c(unique(edge_data$to))
      ))
    
    graph_data <-
      tidygraph::tbl_graph(nodes = node_data,
                           edges = edge_data,
                           directed = FALSE) %>%
      tidygraph::activate(what = "nodes") %>%
      dplyr::mutate(degree = tidygraph::centrality_degree())
    
    g <-
      graph_data %>%
      tidygraph::activate(what = "nodes") %>%
      # dplyr::filter(degree > 3) %>%
      tidygraph::mutate(angle = -360 * (seq_along(name) - 1) / n() + 90)
    
    igraph::V(g)$type <- igraph::bipartite_mapping(g)$type
    
    coords <-
      ggraph::create_layout(g, layout = "bipartite")
    
    coords$index = 1:nrow(coords)
    
    coords$x <-
      coords$x + 1
    
    # coords$y[coords$class == "metabolite"] <- 3 * coords$degree[coords$class == "metabolite"]/max(coords$degree[coords$class == "metabolite"]) + 1
    coords$y[coords$class == "metabolite"] <- 3
    coords$y[coords$class != "metabolite"] <- 1
    
    temp_coords <-
      coords %>%
      dplyr::filter(class != "metabolite") %>%
      dplyr::arrange(x)
    
    temp_coords$x <-
      seq(min(coords$x) + 0.1 * max(coords$x),
          0.9 * max(coords$x),
          length.out = nrow(temp_coords))
    
    coords[coords$class != "metabolite", ] <-
      temp_coords
    
    coords <-
      coords %>%
      dplyr::arrange(index)
    
    coords <-
      coords %>%
      dplyr::mutate(x1 = y, y1 = x) %>%
      dplyr::select(-c(x, y)) %>%
      dplyr::mutate(x = x1, y = y1)
    
    if (circle) {
      coords <-
        coords %>%
        dplyr::mutate(
          theta = y / (max(y) + 1) * 2 * pi,
          r = x + 1,
          x = r * cos(theta),
          y = r * sin(theta)
        )
    }
    
    my_graph <-
      ggraph::create_layout(
        graph = g,
        layout = "manual",
        x = coords$x,
        y = coords$y
        # node.position = coords
      )
    
    plot <-
      ggraph::ggraph(my_graph, layout = 'bipartite') +
      ggraph::geom_edge_diagonal(
        strength = 1,
        aes(color = edge_class),
        edge_width = 0.5,
        alpha = 0.5,
        show.legend = FALSE
      ) +
      ggraph::scale_edge_color_manual(values = c(metabolite_source_color)) +
      ggraph::geom_node_point(aes(size = degree, fill = class),
                              color = "black",
                              shape = 21) +
      ggraph::geom_node_text(
        aes(
          x = x,
          y = y,
          label = ifelse(!class %in% "metabolite", Compound.name, NA)
        ),
        hjust = 1,
        angle = 45,
        size = 4,
        show.legend = FALSE
      ) +
      scale_fill_manual(values = c(metabolite_source_color, metabolite = "black")) +
      ggraph::theme_graph() +
      scale_size_continuous(range = c(2, 10)) +
      theme(plot.background = element_blank(),
            panel.background = element_blank())
    
    plot
    
  }


#' Validate object for metabolite origin analysis
#'
#' @description
#' Checks if an object has the necessary structure and columns for metabolite origin analysis.
#' Validates that required columns for tracking metabolite origins are present.
#'
#' @param object A data.frame or mass_dataset object to validate
#'
#' @return No return value, called for side effects (will stop with error if validation fails)
#'
#' @author Xiaotao Shen \email{xiaotao.shen@@outlook.com}
#'
#' @export

check_object4metablite_origin <-
  function(object) {
    if (!is(object, "data.frame") &
        !is(object, "mass_dataset")) {
      stop("The object should be a data.frame or mass_dataset")
    }
    
    
    if (is(object, "mass_dataset")) {
      object <-
        object@annotation_table
    }
    
    ##check if the object has a column named "HMDB.ID"
    if (all(colnames(object) != "Lab.ID")) {
      stop("The object should have a column named 'Lab.ID'")
    }
    
    ###check the if the object has column named "from_human",
    ###"from_which_part", "from_bacteria", "from_which_bacteria",
    ###"from_fungi", "from_which_fungi", "from_archaea",
    ### "from_which_archaea", "from_plant", "from_which_plant", "from_animal", "from_which_animal"
    ####"from_environment", "from_which_environment" "from_virus", "from_which_virus", "from_protist"
    ###"from_which_protist", "from_drug", "from_which_drug", "from_food", "from_which_food"
    
    if (all(
      colnames(object) %in% c(
        "from_human",
        "from_which_part",
        "from_bacteria",
        "from_which_bacteria",
        "from_fungi",
        "from_which_fungi",
        "from_archaea",
        "from_which_archaea",
        "from_plant",
        "from_which_plant",
        "from_animal",
        "from_which_animal",
        "from_environment",
        "from_which_environment",
        "from_virus",
        "from_which_virus",
        "from_protist",
        "from_which_protist",
        "from_drug",
        "from_which_drug",
        "from_food",
        "from_which_food"
      )
    )) {
      stop(
        "The object should have columns named 'from_human', 'from_which_part', 'from_bacteria', 'from_which_bacteria',
      'from_fungi', 'from_which_fungi', 'from_archaea', 'from_which_archaea', 'from_plant', 'from_which_plant', 'from_animal', 'from_which_animal',
      'from_environment', 'from_which_environment', 'from_virus', 'from_which_virus', 'from_protist', 'from_which_protist', 'from_drug', 'from_which_drug', 'from_food', 'from_which_food'"
      )
    }
    
    
  }


metabolite_source_color <- c(
  "Human" = "#2c61a1",
  "Plant" = "#78938a",
  "Food" = "#f5eddc",
  "Bacteria" = "#0f1c5c",
  "Animal" = "#d2b48c",
  "Enviornment" = "#8f354b",
  "Drug" = "#000000"
)


edge_color <-
  c(
    "metabolite_specific_source" = "#5a435c",
    "source_specific_source" = "#e59589"
  )



#' Create an upset plot of metabolite origins
#'
#' @description
#' Generates an upset plot showing the intersections of metabolites from different sources
#' (human, bacteria, plant, animal, etc.).
#'
#' @param object A mass_dataset object containing metabolite annotation data
#' @param min_size Integer specifying the minimum size of a set to be included in the plot (default: 1)
#' @param counts Logical specifying whether to show counts in the plot (default: TRUE)
#'
#' @return A ggplot2 object representing the upset plot
#'
#' @author Xiaotao Shen \email{xiaotao.shen@@outlook.com}
#'
#' @importFrom ggplot2 ggplot aes geom_bar scale_y_continuous labs theme_bw
#' @importFrom tibble tibble
#' @importFrom purrr map
#' @importFrom dplyr select filter
#' @importFrom stringr str_split
#'
#'
#' @export
metabolite_origin_upsetplot <-
  function(object,
           min_size = 1,
           counts = TRUE) {
    check_object4metablite_origin(object)
    temp_data <-
      object@annotation_table %>%
      dplyr::select(
        Lab.ID,
        from_human,
        from_which_part,
        from_bacteria,
        from_which_bacteria,
        from_plant,
        from_which_plant,
        from_animal,
        from_which_animal,
        from_environment,
        from_which_environment,
        from_drug,
        from_which_drug,
        from_food,
        from_which_food
      )
    
    colnames(temp_data) <-
      c(
        "Lab.ID",
        "Human",
        "Human_name",
        "Bacteria",
        "Bacteria_name",
        "Plant",
        "Plant_name",
        "Animal",
        "Animal_name",
        "Environment",
        "Environment_name",
        "Drug",
        "Drug_name",
        "Food",
        "Food_name"
      )
    
    final_name <-
      c("Human",
        "Bacteria",
        "Plant",
        "Animal",
        "Environment",
        "Drug",
        "Food")
    
    temp_data2 <-
      temp_data[, c(final_name)]
    
    temp_data2[temp_data2 == "Yes"] <- 1
    temp_data2[temp_data2 == "No"] <- 0
    temp_data2[temp_data2 == "Unknown"] <- NA
    temp_data2 <-
      apply(temp_data2, 2, as.numeric) %>%
      as.data.frame()
    
    if (requireNamespace("ComplexUpset", quietly = TRUE)) {
      plot <-
        ComplexUpset::upset(
          data = temp_data2,
          intersect = final_name,
          name = "",
          themes = ComplexUpset::upset_themes,
          width_ratio = 0.15,
          min_size = min_size,
          base_annotations = list('Intersection size' =
                                    ComplexUpset::intersection_size(counts =
                                                                      counts))
        )
    } else {
      stop("Please install the ComplexUpset package")
    }
    
    plot
    
  }



#' Create a network visualization for specific metabolites and their origins
#'
#' @description
#' Generates a detailed network visualization for selected metabolites, showing their
#' connections to source categories and specific sources within those categories.
#'
#' @param object A mass_dataset object containing metabolite annotation data
#' @param metabolite_id Character vector of metabolite IDs to include in the network
#' @param top_specific_source Integer specifying the maximum number of specific sources to show per source category (default: 5)
#'
#' @return A ggplot2 object representing the network visualization
#'
#' @author Xiaotao Shen \email{xiaotao.shen@@outlook.com}
#'
#' @importFrom ggraph ggraph geom_edge_diagonal geom_node_point geom_node_text scale_edge_color_manual theme_graph create_layout
#' @importFrom tidygraph tbl_graph activate mutate
#' @importFrom igraph bipartite_mapping V
#' @importFrom dplyr select filter mutate arrange distinct rename pull left_join all_of
#' @importFrom stringr str_split
#' @importFrom purrr map
#' @importFrom ggplot2 scale_fill_manual theme element_blank
#' @importFrom scales scale_size_continuous
#'
#' @export

metabolite_origin_network <-
  function(object,
           metabolite_id,
           top_specific_source = 5) {
    check_object4metablite_origin(object)
    if (missing(metabolite_id)) {
      stop("Please provide a or several metabolite_ids")
    }
    
    temp_data <-
      object@annotation_table %>%
      dplyr::select(
        Lab.ID,
        from_human,
        from_which_part,
        from_bacteria,
        from_which_bacteria,
        from_plant,
        from_which_plant,
        from_animal,
        from_which_animal,
        from_environment,
        from_which_environment,
        from_drug,
        from_which_drug,
        from_food,
        from_which_food
      ) %>%
      dplyr::filter(Lab.ID %in% metabolite_id)
    
    colnames(temp_data) <-
      c(
        "Lab.ID",
        "Human",
        "Human_name",
        "Bacteria",
        "Bacteria_name",
        "Plant",
        "Plant_name",
        "Animal",
        "Animal_name",
        "Environment",
        "Environment_name",
        "Drug",
        "Drug_name",
        "Food",
        "Food_name"
      )
    
    temp_data$Human_name[temp_data$Human_name == "Unknown"] <-
      paste("Human", "Unknown", sep = "_")
    
    temp_data$Bacteria_name[temp_data$Bacteria_name == "Unknown"] <-
      paste("Bacteria", "Unknown", sep = "_")
    
    temp_data$Plant_name[temp_data$Plant_name == "Unknown"] <-
      paste("Plant", "Unknown", sep = "_")
    
    temp_data$Animal_name[temp_data$Animal_name == "Unknown"] <-
      paste("Animal", "Unknown", sep = "_")
    
    temp_data$Environment_name[temp_data$Environment_name == "Unknown"] <-
      paste("Environment", "Unknown", sep = "_")
    
    temp_data$Drug_name[temp_data$Drug_name == "Unknown"] <-
      paste("Drug", "Unknown", sep = "_")
    
    temp_data$Food_name[temp_data$Food_name == "Unknown"] <-
      paste("Food", "Unknown", sep = "_")
    
    Source <-
      unique(temp_data$Lab.ID) %>%
      purrr::map(function(x) {
        temp <-
          temp_data %>%
          dplyr::filter(Lab.ID == x)
        colnames(temp)[which(as.character(temp) == "Yes")]
      })
    
    temp_data2 <-
      temp_data %>%
      dplyr::select(Lab.ID)
    
    temp_data2$source <- Source
    
    ###source information
    temp_data3 <-
      seq_len(nrow(temp_data2)) %>%
      purrr::map(function(i) {
        data.frame(Lab.ID = temp_data2$Lab.ID[i], Source = temp_data2$source[[i]])
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    #####specific source
    temp_data4 <-
      unique(temp_data3$Lab.ID) %>%
      purrr::map(function(x) {
        temp <-
          temp_data3 %>%
          dplyr::filter(Lab.ID == x)
        index <-
          temp %>%
          dplyr::pull(Source) %>%
          paste0("_name")
        
        specific_source <-
          temp_data %>%
          dplyr::filter(Lab.ID == x) %>%
          dplyr::select(all_of(index)) %>%
          as.character() %>%
          purrr::map(function(y) {
            temp <-
              stringr::str_split(y, "\\{\\}")[[1]]
            head(temp, top_specific_source)
          })
        
        temp <-
          seq_len(nrow(temp)) %>%
          purrr::map(function(i) {
            data.frame(
              Lab.ID = temp$Lab.ID[i],
              Source = temp$Source[i],
              Specific_source = specific_source[[i]]
            )
          }) %>%
          do.call(rbind, .) %>%
          as.data.frame()
        
        temp
        
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    
    edge_data1 <-
      temp_data4 %>%
      dplyr::select(Lab.ID, Specific_source) %>%
      dplyr::rename(from = Lab.ID, to = Specific_source) %>%
      dplyr::distinct(from, to, .keep_all = TRUE) %>%
      dplyr::mutate(edge_class = "metabolite_specific_source")
    
    edge_data2 <-
      temp_data4 %>%
      dplyr::select(Source, Specific_source) %>%
      dplyr::rename(from = Source, to = Specific_source) %>%
      dplyr::distinct(from, to, .keep_all = TRUE) %>%
      dplyr::mutate(edge_class = from)
    
    edge_data <-
      rbind(edge_data1, edge_data2)
    
    node_data <-
      rbind(
        data.frame(
          name = temp_data4$Lab.ID,
          node_class = "metabolite",
          node_class2 = "metabolite"
        ),
        data.frame(
          name = temp_data4$Source,
          node_class = "source",
          node_class2 = temp_data4$Source
        ),
        data.frame(
          name = temp_data4$Specific_source,
          node_class = "specific_source",
          node_class2 = temp_data4$Source
        )
      ) %>%
      dplyr::distinct(name, node_class, .keep_all = TRUE)
    
    
    node_data <-
      node_data %>%
      dplyr::left_join(object@annotation_table[, c("Lab.ID", "Compound.name")], by = c("name" = "Lab.ID")) %>%
      dplyr::mutate(node_name = ifelse(node_class == "metabolite", Compound.name, name))
    
    graph_data <-
      tidygraph::tbl_graph(nodes = node_data,
                           edges = edge_data,
                           directed = FALSE) %>%
      tidygraph::activate(what = "nodes") %>%
      dplyr::mutate(degree = tidygraph::centrality_degree())
    
    g <-
      graph_data %>%
      tidygraph::activate(what = "nodes") %>%
      # dplyr::filter(degree > 3) %>%
      tidygraph::mutate(angle = -360 * (seq_along(name) - 1) / n() + 90)
    
    igraph::V(g)$type <- igraph::bipartite_mapping(g)$type
    
    coords <-
      ggraph::create_layout(g, layout = "bipartite")
    
    coords$index = 1:nrow(coords)
    
    coords$x <-
      coords$x + 1
    
    coords$y[coords$node_class == "metabolite"] <- 3
    coords$y[coords$node_class == "specific_source"] <- 2
    coords$y[coords$node_class == "source"] <- 1
    
    temp_coords <-
      coords %>%
      dplyr::filter(node_class == "source") %>%
      dplyr::arrange(x)
    
    temp_coords$x <-
      seq(min(coords$x) + 0.1 * max(coords$x),
          0.9 * max(coords$x),
          length.out = nrow(temp_coords))
    
    coords[coords$node_class == "source", ] <-
      temp_coords
    
    coords <-
      coords %>%
      dplyr::arrange(index)
    
    coords <-
      coords %>%
      dplyr::mutate(x1 = y, y1 = x) %>%
      dplyr::select(-c(x, y)) %>%
      dplyr::mutate(x = x1, y = y1)
    
    my_graph <-
      ggraph::create_layout(
        graph = g,
        layout = "manual",
        x = coords$x,
        y = coords$y
      )
    
    plot <-
      ggraph::ggraph(my_graph, layout = 'bipartite') +
      ggraph::geom_edge_diagonal(
        strength = 1,
        aes(color = edge_class),
        edge_width = 1,
        alpha = 0.7,
        show.legend = FALSE
      ) +
      ggraph::scale_edge_color_manual(values = c(
        metabolite_source_color,
        metabolite_specific_source = "#5a435c"
      )) +
      ggraph::geom_node_point(aes(size = degree, fill = node_class2),
                              color = "black",
                              shape = 21) +
      ggraph::geom_node_text(
        aes(x = x, y = y, label = node_name),
        hjust = 1,
        angle = 45,
        size = 4,
        show.legend = FALSE
      ) +
      scale_fill_manual(
        values = c(
          metabolite_source_color,
          metabolite = "black",
          specific_source = "#2dbc94"
        )
      ) +
      ggraph::theme_graph() +
      scale_size_continuous(range = c(2, 10)) +
      theme(plot.background = element_blank(),
            panel.background = element_blank())
    
    plot
    
  }


#' Create a network visualization for metabolites from a specific source
#'
#' @description
#' Generates a network visualization showing metabolites originating from a specified
#' source category (e.g., Human, Bacteria, Plant).
#'
#' @param object A mass_dataset object containing metabolite annotation data
#' @param source_id Character vector of source categories to include (e.g., "Human", "Bacteria", "Plant")
#' @param top_specific_source Integer specifying the maximum number of specific sources to show per source category (default: 5)
#'
#' @return A ggplot2 object representing the network visualization
#'
#' @author Xiaotao Shen \email{xiaotao.shen@@outlook.com}
#'
#' @importFrom ggraph ggraph geom_edge_diagonal geom_node_point geom_node_text scale_edge_color_manual theme_graph create_layout
#' @importFrom tidygraph tbl_graph activate mutate
#' @importFrom igraph bipartite_mapping V
#' @importFrom dplyr select filter mutate arrange distinct rename pull left_join all_of
#' @importFrom stringr str_split
#' @importFrom purrr map
#' @importFrom ggplot2 scale_fill_manual theme element_blank
#' @importFrom scales scale_size_continuous
#'
#' @export

source_network <-
  function(object, source_id, top_specific_source = 5) {
    check_object4metablite_origin(object)
    if (missing(source_id)) {
      stop(
        "Please provide a or several source_ids, such as Human, Bacteria, Plant, Animal, Environment, Drug, Food"
      )
    }
    
    temp_data <-
      object@annotation_table %>%
      dplyr::select(
        Lab.ID,
        from_human,
        from_which_part,
        from_bacteria,
        from_which_bacteria,
        from_plant,
        from_which_plant,
        from_animal,
        from_which_animal,
        from_environment,
        from_which_environment,
        from_drug,
        from_which_drug,
        from_food,
        from_which_food
      )
    
    colnames(temp_data) <-
      c(
        "Lab.ID",
        "Human",
        "Human_name",
        "Bacteria",
        "Bacteria_name",
        "Plant",
        "Plant_name",
        "Animal",
        "Animal_name",
        "Environment",
        "Environment_name",
        "Drug",
        "Drug_name",
        "Food",
        "Food_name"
      )
    
    temp_data$Human_name[temp_data$Human_name == "Unknown"] <-
      paste("Human", "Unknown", sep = "_")
    
    temp_data$Bacteria_name[temp_data$Bacteria_name == "Unknown"] <-
      paste("Bacteria", "Unknown", sep = "_")
    
    temp_data$Plant_name[temp_data$Plant_name == "Unknown"] <-
      paste("Plant", "Unknown", sep = "_")
    
    temp_data$Animal_name[temp_data$Animal_name == "Unknown"] <-
      paste("Animal", "Unknown", sep = "_")
    
    temp_data$Environment_name[temp_data$Environment_name == "Unknown"] <-
      paste("Environment", "Unknown", sep = "_")
    
    temp_data$Drug_name[temp_data$Drug_name == "Unknown"] <-
      paste("Drug", "Unknown", sep = "_")
    
    temp_data$Food_name[temp_data$Food_name == "Unknown"] <-
      paste("Food", "Unknown", sep = "_")
    
    Source <-
      unique(temp_data$Lab.ID) %>%
      purrr::map(function(x) {
        temp <-
          temp_data %>%
          dplyr::filter(Lab.ID == x)
        colnames(temp)[which(as.character(temp) == "Yes")]
      })
    
    temp_data2 <-
      temp_data %>%
      dplyr::select(Lab.ID)
    
    temp_data2$source <- Source
    
    ###source information
    temp_data3 <-
      seq_len(nrow(temp_data2)) %>%
      purrr::map(function(i) {
        temp <-
          tryCatch(
            data.frame(Lab.ID = temp_data2$Lab.ID[i], Source = temp_data2$source[[i]]),
            error = function(e) {
              NULL
            }
          )
        temp
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame() %>%
      dplyr::filter(Source %in% source_id)
    
    #####specific source
    temp_data4 <-
      unique(temp_data3$Lab.ID) %>%
      purrr::map(function(x) {
        temp <-
          temp_data3 %>%
          dplyr::filter(Lab.ID == x)
        index <-
          temp %>%
          dplyr::pull(Source) %>%
          paste0("_name")
        
        specific_source <-
          temp_data %>%
          dplyr::filter(Lab.ID == x) %>%
          dplyr::select(all_of(index)) %>%
          as.character() %>%
          purrr::map(function(y) {
            temp <-
              stringr::str_split(y, "\\{\\}")[[1]]
            head(temp, top_specific_source)
          })
        
        temp <-
          seq_len(nrow(temp)) %>%
          purrr::map(function(i) {
            data.frame(
              Lab.ID = temp$Lab.ID[i],
              Source = temp$Source[i],
              Specific_source = specific_source[[i]]
            )
          }) %>%
          do.call(rbind, .) %>%
          as.data.frame()
        
        temp
        
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame() %>%
      dplyr::filter(Source %in% source_id)
    
    
    edge_data1 <-
      temp_data4 %>%
      dplyr::select(Lab.ID, Specific_source) %>%
      dplyr::rename(from = Lab.ID, to = Specific_source) %>%
      dplyr::distinct(from, to, .keep_all = TRUE) %>%
      dplyr::mutate(edge_class = "metabolite_specific_source")
    
    edge_data2 <-
      temp_data4 %>%
      dplyr::select(Source, Specific_source) %>%
      dplyr::rename(from = Source, to = Specific_source) %>%
      dplyr::distinct(from, to, .keep_all = TRUE) %>%
      dplyr::mutate(edge_class = from)
    
    edge_data <-
      rbind(edge_data1, edge_data2)
    
    node_data <-
      rbind(
        data.frame(
          name = temp_data4$Lab.ID,
          node_class = "metabolite",
          node_class2 = "metabolite"
        ),
        data.frame(
          name = temp_data4$Source,
          node_class = "source",
          node_class2 = temp_data4$Source
        ),
        data.frame(
          name = temp_data4$Specific_source,
          node_class = "specific_source",
          node_class2 = temp_data4$Source
        )
      ) %>%
      dplyr::distinct(name, node_class, .keep_all = TRUE)
    
    
    node_data <-
      node_data %>%
      dplyr::left_join(object@annotation_table[, c("Lab.ID", "Compound.name")], by = c("name" = "Lab.ID")) %>%
      dplyr::mutate(node_name = ifelse(node_class == "metabolite", Compound.name, name))
    
    graph_data <-
      tidygraph::tbl_graph(nodes = node_data,
                           edges = edge_data,
                           directed = FALSE) %>%
      tidygraph::activate(what = "nodes") %>%
      dplyr::mutate(degree = tidygraph::centrality_degree())
    
    g <-
      graph_data %>%
      tidygraph::activate(what = "nodes") %>%
      # dplyr::filter(degree > 3) %>%
      tidygraph::mutate(angle = -360 * (seq_along(name) - 1) / n() + 90)
    
    igraph::V(g)$type <- igraph::bipartite_mapping(g)$type
    
    coords <-
      ggraph::create_layout(g, layout = "bipartite")
    
    coords$index = 1:nrow(coords)
    
    coords$x <-
      coords$x + 1
    
    coords$y[coords$node_class == "metabolite"] <- 3
    coords$y[coords$node_class == "specific_source"] <- 2
    coords$y[coords$node_class == "source"] <- 1
    
    temp_coords <-
      coords %>%
      dplyr::filter(node_class == "source") %>%
      dplyr::arrange(x)
    
    temp_coords$x <-
      seq(min(coords$x) + 0.1 * max(coords$x),
          0.9 * max(coords$x),
          length.out = nrow(temp_coords))
    
    coords[coords$node_class == "source", ] <-
      temp_coords
    
    coords <-
      coords %>%
      dplyr::arrange(index)
    
    coords <-
      coords %>%
      dplyr::mutate(x1 = y, y1 = x) %>%
      dplyr::select(-c(x, y)) %>%
      dplyr::mutate(x = x1, y = y1)
    
    my_graph <-
      ggraph::create_layout(
        graph = g,
        layout = "manual",
        x = coords$x,
        y = coords$y
      )
    
    plot <-
      ggraph::ggraph(my_graph, layout = 'bipartite') +
      ggraph::geom_edge_diagonal(
        strength = 1,
        aes(color = edge_class),
        edge_width = 1,
        alpha = 0.7,
        show.legend = FALSE
      ) +
      ggraph::scale_edge_color_manual(values = c(
        metabolite_source_color,
        metabolite_specific_source = "#5a435c"
      )) +
      ggraph::geom_node_point(aes(size = degree, fill = node_class2),
                              color = "black",
                              shape = 21) +
      ggraph::geom_node_text(
        aes(x = x, y = y, label = node_name),
        hjust = 1,
        angle = 45,
        size = 4,
        show.legend = FALSE
      ) +
      scale_fill_manual(
        values = c(
          metabolite_source_color,
          metabolite = "black",
          specific_source = "#2dbc94"
        )
      ) +
      ggraph::theme_graph() +
      scale_size_continuous(range = c(2, 10)) +
      theme(plot.background = element_blank(),
            panel.background = element_blank())
    
    plot
    
  }


#' Create a network visualization for metabolites from a specific sub-source
#'
#' @description
#' Generates a network visualization showing metabolites originating from a specified
#' specific source (e.g., a specific tissue, bacterial species, etc.).
#'
#' @param object A mass_dataset object containing metabolite annotation data
#' @param specific_source_id Character vector of specific source identifiers to include
#' @param top_specific_source Integer specifying the maximum number of specific sources to show per source category (default: 5)
#'
#' @return A ggplot2 object representing the network visualization
#'
#' @author Xiaotao Shen \email{xiaotao.shen@@outlook.com}
#'
#' @importFrom ggraph ggraph geom_edge_diagonal geom_node_point geom_node_text scale_edge_color_manual theme_graph create_layout
#' @importFrom tidygraph tbl_graph activate mutate
#' @importFrom igraph bipartite_mapping V
#' @importFrom dplyr select filter mutate arrange distinct rename pull left_join all_of
#' @importFrom stringr str_split
#' @importFrom purrr map
#' @importFrom ggplot2 scale_fill_manual theme element_blank
#' @importFrom scales scale_size_continuous
#'
#' @export
specific_source_network <-
  function(object,
           specific_source_id,
           top_specific_source = 5) {
    check_object4metablite_origin(object)
    if (missing(specific_source_id)) {
      stop("Please provide a or several specific_source_id")
    }
    
    temp_data <-
      object@annotation_table %>%
      dplyr::select(
        Lab.ID,
        from_human,
        from_which_part,
        from_bacteria,
        from_which_bacteria,
        from_plant,
        from_which_plant,
        from_animal,
        from_which_animal,
        from_environment,
        from_which_environment,
        from_drug,
        from_which_drug,
        from_food,
        from_which_food
      )
    
    colnames(temp_data) <-
      c(
        "Lab.ID",
        "Human",
        "Human_name",
        "Bacteria",
        "Bacteria_name",
        "Plant",
        "Plant_name",
        "Animal",
        "Animal_name",
        "Environment",
        "Environment_name",
        "Drug",
        "Drug_name",
        "Food",
        "Food_name"
      )
    
    temp_data$Human_name[temp_data$Human_name == "Unknown"] <-
      paste("Human", "Unknown", sep = "_")
    
    temp_data$Bacteria_name[temp_data$Bacteria_name == "Unknown"] <-
      paste("Bacteria", "Unknown", sep = "_")
    
    temp_data$Plant_name[temp_data$Plant_name == "Unknown"] <-
      paste("Plant", "Unknown", sep = "_")
    
    temp_data$Animal_name[temp_data$Animal_name == "Unknown"] <-
      paste("Animal", "Unknown", sep = "_")
    
    temp_data$Environment_name[temp_data$Environment_name == "Unknown"] <-
      paste("Environment", "Unknown", sep = "_")
    
    temp_data$Drug_name[temp_data$Drug_name == "Unknown"] <-
      paste("Drug", "Unknown", sep = "_")
    
    temp_data$Food_name[temp_data$Food_name == "Unknown"] <-
      paste("Food", "Unknown", sep = "_")
    
    Source <-
      unique(temp_data$Lab.ID) %>%
      purrr::map(function(x) {
        temp <-
          temp_data %>%
          dplyr::filter(Lab.ID == x)
        colnames(temp)[which(as.character(temp) == "Yes")]
      })
    
    temp_data2 <-
      temp_data %>%
      dplyr::select(Lab.ID)
    
    temp_data2$source <- Source
    
    ###source information
    temp_data3 <-
      seq_len(nrow(temp_data2)) %>%
      purrr::map(function(i) {
        temp <-
          tryCatch(
            data.frame(Lab.ID = temp_data2$Lab.ID[i], Source = temp_data2$source[[i]]),
            error = function(e) {
              NULL
            }
          )
        temp
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    #####specific source
    temp_data4 <-
      unique(temp_data3$Lab.ID) %>%
      purrr::map(function(x) {
        temp <-
          temp_data3 %>%
          dplyr::filter(Lab.ID == x)
        index <-
          temp %>%
          dplyr::pull(Source) %>%
          paste0("_name")
        
        specific_source <-
          temp_data %>%
          dplyr::filter(Lab.ID == x) %>%
          dplyr::select(all_of(index)) %>%
          as.character() %>%
          purrr::map(function(y) {
            temp <-
              stringr::str_split(y, "\\{\\}")[[1]]
            head(temp, top_specific_source)
          })
        
        temp <-
          seq_len(nrow(temp)) %>%
          purrr::map(function(i) {
            data.frame(
              Lab.ID = temp$Lab.ID[i],
              Source = temp$Source[i],
              Specific_source = specific_source[[i]]
            )
          }) %>%
          do.call(rbind, .) %>%
          as.data.frame()
        
        temp
        
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame() %>%
      dplyr::filter(Specific_source %in% specific_source_id)
    
    edge_data1 <-
      temp_data4 %>%
      dplyr::select(Lab.ID, Specific_source) %>%
      dplyr::rename(from = Lab.ID, to = Specific_source) %>%
      dplyr::distinct(from, to, .keep_all = TRUE) %>%
      dplyr::mutate(edge_class = "metabolite_specific_source")
    
    edge_data2 <-
      temp_data4 %>%
      dplyr::select(Source, Specific_source) %>%
      dplyr::rename(from = Source, to = Specific_source) %>%
      dplyr::distinct(from, to, .keep_all = TRUE) %>%
      dplyr::mutate(edge_class = from)
    
    edge_data <-
      rbind(edge_data1, edge_data2)
    
    node_data <-
      rbind(
        data.frame(
          name = temp_data4$Lab.ID,
          node_class = "metabolite",
          node_class2 = "metabolite"
        ),
        data.frame(
          name = temp_data4$Source,
          node_class = "source",
          node_class2 = temp_data4$Source
        ),
        data.frame(
          name = temp_data4$Specific_source,
          node_class = "specific_source",
          node_class2 = temp_data4$Source
        )
      ) %>%
      dplyr::distinct(name, node_class, .keep_all = TRUE)
    
    
    node_data <-
      node_data %>%
      dplyr::left_join(object@annotation_table[, c("Lab.ID", "Compound.name")], by = c("name" = "Lab.ID")) %>%
      dplyr::mutate(node_name = ifelse(node_class == "metabolite", Compound.name, name))
    
    graph_data <-
      tidygraph::tbl_graph(nodes = node_data,
                           edges = edge_data,
                           directed = FALSE) %>%
      tidygraph::activate(what = "nodes") %>%
      dplyr::mutate(degree = tidygraph::centrality_degree())
    
    g <-
      graph_data %>%
      tidygraph::activate(what = "nodes") %>%
      # dplyr::filter(degree > 3) %>%
      tidygraph::mutate(angle = -360 * (seq_along(name) - 1) / n() + 90)
    
    igraph::V(g)$type <- igraph::bipartite_mapping(g)$type
    
    coords <-
      ggraph::create_layout(g, layout = "bipartite")
    
    coords$index = 1:nrow(coords)
    
    coords$x <-
      coords$x + 1
    
    coords$y[coords$node_class == "metabolite"] <- 3
    coords$y[coords$node_class == "specific_source"] <- 2
    coords$y[coords$node_class == "source"] <- 1
    
    temp_coords <-
      coords %>%
      dplyr::filter(node_class == "source") %>%
      dplyr::arrange(x)
    
    temp_coords$x <-
      seq(min(coords$x) + 0.1 * max(coords$x),
          0.9 * max(coords$x),
          length.out = nrow(temp_coords))
    
    coords[coords$node_class == "source", ] <-
      temp_coords
    
    coords <-
      coords %>%
      dplyr::arrange(index)
    
    coords <-
      coords %>%
      dplyr::mutate(x1 = y, y1 = x) %>%
      dplyr::select(-c(x, y)) %>%
      dplyr::mutate(x = x1, y = y1)
    
    my_graph <-
      ggraph::create_layout(
        graph = g,
        layout = "manual",
        x = coords$x,
        y = coords$y
      )
    
    plot <-
      ggraph::ggraph(my_graph, layout = 'bipartite') +
      ggraph::geom_edge_diagonal(
        strength = 1,
        aes(color = edge_class),
        edge_width = 1,
        alpha = 0.7,
        show.legend = FALSE
      ) +
      ggraph::scale_edge_color_manual(values = c(
        metabolite_source_color,
        metabolite_specific_source = "#5a435c"
      )) +
      ggraph::geom_node_point(aes(size = degree, fill = node_class2),
                              color = "black",
                              shape = 21) +
      ggraph::geom_node_text(
        aes(x = x, y = y, label = node_name),
        hjust = 1,
        angle = 45,
        size = 4,
        show.legend = FALSE
      ) +
      scale_fill_manual(
        values = c(
          metabolite_source_color,
          metabolite = "black",
          specific_source = "#2dbc94"
        )
      ) +
      ggraph::theme_graph() +
      scale_size_continuous(range = c(2, 10)) +
      theme(plot.background = element_blank(),
            panel.background = element_blank())
    
    plot
    
  }
