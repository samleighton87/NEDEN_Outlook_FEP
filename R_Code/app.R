#load required libraries
library(shiny)
library(PostcodesioR)
library(caret)
library(png)
library(readr)
library(RANN)

#load objects and data
ccg_local_conc_2019 = read_csv("CCG_Local_Conc_2019.csv")
load("mod_eden_all_final_Rem_outlook_glm_new_shrink.rda")

#Define custom methods
smileyDiagram = function(probRaw, estimated, postCode, CCG)
{
  probPercent = probRaw*100
  if(probPercent > 100 | probPercent < 0)
  {
    print("Enter a percentage probability between 0 and 100")
  }else
  {
    #Smiley face diagram
    #data with fixed ratio and random distribution
    probPercent = floor(probPercent + 0.5)
    x=sample(as.factor(rep(c("A","B"),c(as.integer(probPercent), as.integer(100-probPercent)))))
    # input parameters - nr * nc should equal length(x)
    
    #colours are actually only used for mapping - they just have to be called a valid colour
    #here blue is happy and red is sad but the colours only change if you change the pictures
    cols = NULL
    if(probPercent == 0)
    {
      cols <- c("red", "blue")
    }else
    {
      cols <- c("blue", "red")
    }
    
    nr <- 10
    nc <- 10
    # create data.frame of positions and colors
    m <- matrix(cols[x], nr, nc)
    DF <- data.frame(row = c(row(m)), col = c(col(m)[, nc:1]), value = c(m), 
                     stringsAsFactors = FALSE)
    
    # blank graph to insert man icons into
    xp <- 1.25
    if(!estimated)
    {
      plot(col ~ row, DF, col = DF$value, asp = 1,
           xlim = c(0, xp * nr), ylim = c(0, xp * nc),
           axes = FALSE, xlab = "", ylab = "", type = "n",
           main = paste0(floor((probRaw*1000) + 0.5)/10, "% risk of not being in remission at 1 year\n\n"),
           sub = paste0("In a crowd of 100 people with the same\nrisk factors, ", probPercent, " are likely to not be\nin remission from FEP at 1 year\n\n(Deprivation score for ",CCG,")"),
           cex.main=1.5, cex.sub =1.25)
    }else if(!is.na(CCG))
    {
      plot(col ~ row, DF, col = DF$value, asp = 1,
           xlim = c(0, xp * nr), ylim = c(0, xp * nc),
           axes = FALSE, xlab = "", ylab = "", type = "n",
           main = paste0(floor((probRaw*1000) + 0.5)/10, "% risk of not being in remission at 1 year.\n(Some data estimated as values left blank)\n"),
           sub = paste0("In a crowd of 100 people with the same\nrisk factors, ", probPercent, " are likely to not be\nin remission from FEP at 1 year.\n\n(Deprivation score for ",CCG,")"),
                        col.sub = "black",cex.main=1.5, cex.sub =1.25)
    }else
    {
      plot(col ~ row, DF, col = DF$value, asp = 1,
           xlim = c(0, xp * nr), ylim = c(0, xp * nc),
           axes = FALSE, xlab = "", ylab = "", type = "n",
           main = paste0(floor((probRaw*1000) + 0.5)/10, "% risk of not being in remission at 1 year.\n(Some data estimated as values left blank)\n"),
           sub = paste0("In a crowd of 100 people with the same\nrisk factors, ", probPercent, " are likely to not be\nin remission from FEP at 1 year.\n\n(No or invalid postcode - deprivation estimated)"),
           col.sub = "black",cex.main=1.5, cex.sub =1.25)
    }
    
    #images must be in working directory
    smile_face <- readPNG("smile_face.png")
    sad_face <- readPNG("sad_face.png")
    
    G <- subset(transform(DF, row = xp * row, col = xp * col), value == "blue")
    with(G, rasterImage(sad_face,
                        row - .5, col - .5, row + .5, col + .5, 
                        xlim = c(0, xp * nr), ylim = c(0, xp * nc),
                        xlab = "", ylab = ""))
    
    R <- subset(transform(DF, row = xp * row, col = xp * col), value == "red")
    with(R, rasterImage(smile_face,
                        row - .5, col - .5, row + .5, col + .5, 
                        xlim = c(0, xp * nr), ylim = c(0, xp * nc),
                        xlab = "", ylab = ""))
  }
}

#takes postcode and 4 column table of ccg code, ccg name, raw imd score and standardised score
getCCGFromPostcode = function(postcode, ccgTable)
{
  letters_only <- function(x) !grepl("[^A-Za-z]", x)
  numbers_only <- function(x) !grepl("\\D", x)
  outCodeToPostCode = random_postcode(postcode)
  completePostcodeBool = postcode_validation(postcode)
  completePostCode = NULL
  CCGName = NULL
  DepZScore = NULL
  if(!is.null(outCodeToPostCode))
  {
    #is a valid outcode
    print("Valid partial postcode outcode")
    if(outCodeToPostCode$country != "England")
    {
      print("Sorry only valid for NHS England")
    }else
    {
      ccg = ccgTable[which(ccgTable[[1]] == outCodeToPostCode$codes$ccg),]
      print(paste0("Your postocde outcode corresponds to ",ccg[[2]]))
      print(paste0("The Deprivation Z-Score for ",ccg[[2]]," is ", ccg[[4]], "."))
      completePostCode = outCodeToPostCode$postcode
      CCGName = ccg[[2]]
      DepZScore = ccg[[4]]
    }
  }else if (completePostcodeBool)
  {
    print("Valid complete postcode")
    #is a valid complete postcode
    completePostCode = postcode_lookup(postcode)
    if(completePostCode$country != "England")
    {
      print("Sorry only valid for NHS England")
      completePostCode = NULL
    }else
    {
      ccg = ccgTable[which(ccgTable[[1]] ==  completePostCode$ccg_code),]
      print(paste0("Your postocde corresponds to ",ccg[[2]]))
      print(paste0("The Deprivation Z-Score for ",ccg[[2]]," is ", ccg[[4]], "."))
      CCGName = ccg[[2]]
      DepZScore = ccg[[4]]
    }
  }else if (!letters_only(postcode) & #not all letters
            !numbers_only(postcode) & #not all numbers
            #don't autocomplete partial outcodes 
            #so at least 3 characters if second character is a number (if it is a complete outcode it is handled above)
            ((nchar(postcode)>2 & numbers_only(strsplit(as.character(postcode), "")[[1]][2]))|
             #or at least 4 characters if second character is a letter (if it is a complete outcode it is handled above)
             (nchar(postcode)>3 & letters_only(strsplit(as.character(postcode), "")[[1]][2])))
  )
  {
    tryCatch({#not an outcode but partial postcode greater than outcode
      autoCompletePostCode = postcode_autocomplete(postcode, 1)[[1]]
      completePostCode = postcode_lookup(autoCompletePostCode)
      print("Valid partial postcode longer than outcode")
      if(completePostCode$country != "England")
      {
        print("Sorry only valid for NHS England")
        completePostCode = NULL
      }else
      {
        ccg = ccgTable[which(ccgTable[[1]] ==  completePostCode$ccg_code),]
        print(paste0("Your partial postocde corresponds to ",ccg[[2]]))
        print(paste0("The Deprivation Z-Score for ",ccg[[2]]," is ", ccg[[4]], "."))
        CCGName = ccg[[2]]
        DepZScore = ccg[[4]]
      }
    },
    #not a postcode
    error=function(error_message) {
      print("Neither a valid complete postcode nor a valid partial postcode")
    })
  }else
  {
    #too short and meaningless
    print("Please enter at least a valid postcode outcode")
  }
  if(is.null(completePostCode))
  {
    x = list("completePostCode" = NA, "CCGName" = NA, "DepZScore" = NA)
  }else
  {
    x = list("completePostCode" = completePostCode, "CCGName" = CCGName, "DepZScore" = DepZScore)
    return(x)
  }
}

ui <- fluidPage(
  HTML('<meta name="viewport" content="width=1024">'),
  
  titlePanel("Calculate the Risk of Not Being in Remission from a First Episode of Psychosis at 1 Year."),
  sidebarLayout(sidebarPanel(
    h2("Enter Baseline Variables (if known):"),
    br(),
    selectInput("Qualification_Level_Ordinal","Highest Qualification Level",
                choices = list(" " = -1,"None" = 0, "GCSE/NVQ level 1 or 2" = 1,
                               "A-level/GNVQ/BTEC/NVQ level 3" = 2,
                               "Degree/HND/NVQ level 4 or above" = 3), selected = -1),
    numericInput("ADJ_DUP","Duration of Untreated Psychosis (days)", value = NA),
    selectInput("Client_Sociability_Withdrawal_Late_Adolescence","PAS Late Adolescence - Sociability & Withdrawal",
                choices = list(" " = -1,"0 - Not withdrawn" = 0, "1" = 1,"2" = 2,"3" = 3,"4" = 4, "5" = 5, "6 - Withdrawn & isolated"=6), selected = -1),
    selectInput("Client_Highest_Functioning_Achieved_Life","PAS General - Degree of Interest in Life",
                choices = list(" " = -1,"0 - Fully able to function" = 0, "1" = 1,"2" = 2,"3" = 3,"4" = 4, "5" = 5, "6 - Unable to function"=6), selected = -1),
    selectInput("Client_Energy_Level","PAS General - Energy Level",
                choices = list(" " = -1 ,"0 - Strong drive, keen, active, alert" = 0, "1" = 1,"2" = 2,"3" = 3,"4" = 4, "5" = 5, "6 - Submissive, inadequate, passive"=6), selected = -1),
    selectInput("BL_PANSS_P3_Hallucinatory_Behaviour","PANSS P3 - Hallucinatory Behaviour",
                choices = list(" " = -1, "1 - Absent" = 1,"2"= 2,"3" = 3,"4" = 4, "5" = 5,"6"=6, "7 - Extreme"=7), selected = -1),
    selectInput("BL_PANSS_N4_Passive_Social_Withdrawal","PANSS P4 - Passive Social Withdrawal",
                choices = list(" " = -1, "1 - Absent" = 1,"2"= 2,"3" = 3,"4" = 4, "5" = 5,"6"=6, "7 - Extreme"=7), selected = -1),
    selectInput("BL_PANSS_G9_Unusual_Thought_Content","PANSS G9 - Unusual Thought Content",
                choices = list(" " = -1, "1 - Absent" = 1,"2"= 2,"3" = 3,"4" = 4, "5" = 5,"6"=6, "7 - Extreme"=7), selected = -1),
    selectInput("YMRS_Appearance","YMRS - Appearance",
                choices = list(" " = -1,"0 - Appropriate dress & grooming" = 0,"1" = 1,"2" = 2,"3" = 3,"4 - Completely unkempt, decorated, bizarre garb" = 4), selected = -1),
    selectInput("IS_Not_Due_Illness","Insight Scale - None of the unusual things I experienced are due to an illness",
                choices = list(" " = -1,"0 - Agree" = 0,"1 - Unsure" = 1,"2 - Disagree" = 2), selected = -1),
    numericInput("BL_GAF_Symptoms","Global Assessment of Functioning - Symptoms (1 - Persistent to 100 - Absent)", value = NA, min = 1, max = 100),
    textInput("PostCodeToDep","Postcode (enter at least postcode district, spaces are ignored) - England only, otherwise deprivation estimated",value = " ")
  ),
  mainPanel(
    br(),br(),actionButton("goButton", "Calculate!"),br(),br(),
    plotOutput("probNoRemission", width = "100%"),
    br(),br(),actionButton("resetButton", "Reset"),br(),br()
  )
))

# Define server logic required to draw a piechart ----
server <- function(input, output, session) {

  output$probNoRemission = renderPlot({
    plot(0,type='n',axes=FALSE, main="Welcome to the First Episode Psychosis\nNon-Remission Calculator v0.0.2 (25/11/19)", 
    sub="We used data from the National Evaluation of Development of Early 
          intervention Network study (NEDEN) for model development and 
          internal-external validation. NEDEN is a longitudinal naturalistic 
          study of 1027 patients aged 14-35 with first-episode psychosis
          recruited from 14 early intervention services across the National 
          Health Service (NHS) in England (2005-10).

        For external validation, we used data from the Outlook longitudinal 
        naturalistic study of 399 patients aged 16-35 recruited from an
        additional 11 NHS England early intervention services, throughout 
        April 2006-February 2009."
        , ylab = "",xlab = "",cex.main = 2, col.main = "black")
    text(x = 1, y=0.5, "Developed in collaboration between the Institute of\nMental Health & Wellbeing, University of Glasgow &\nInstitute of Mental Health, University of Birmingham\nby Dr Samuel Leighton and Dr Pavan Mallikarjun.", cex = 1.5)
  })
  
  observeEvent(input$BL_GAF_Symptoms, {
    if(as.numeric(input$BL_GAF_Symptoms) < 1 | as.numeric(input$BL_GAF_Symptoms) > 100 | is.na(input$BL_GAF_Symptoms))
    {
      updateNumericInput(session = session, inputId = "BL_GAF_Symptoms",value = NA)
    }else
    {
      return()
    }
  })
  
  observeEvent(input$resetButton, {
      updateSelectInput(session = session, inputId = "Qualification_Level_Ordinal",selected = -1)
      updateNumericInput(session = session, inputId = "ADJ_DUP", value = NA)
      updateSelectInput(session = session, inputId = "Client_Sociability_Withdrawal_Late_Adolescence",selected = -1)
      updateSelectInput(session = session, inputId = "Client_Highest_Functioning_Achieved_Life",selected = -1)
      updateSelectInput(session = session, inputId = "Client_Energy_Level",selected = -1) 
      updateSelectInput(session = session, inputId = "BL_PANSS_P3_Hallucinatory_Behaviour",selected = -1)
      updateSelectInput(session = session, inputId = "BL_PANSS_N4_Passive_Social_Withdrawal",selected = -1)
      updateSelectInput(session = session, inputId = "BL_PANSS_G9_Unusual_Thought_Content",selected = -1)
      updateSelectInput(session = session, inputId = "YMRS_Appearance",selected = -1)
      updateSelectInput(session = session, inputId = "IS_Not_Due_Illness",selected = -1)
      updateNumericInput(session = session, inputId = "BL_GAF_Symptoms",value = NA)
      updateTextInput(session = session, inputId = "PostCodeToDep",value = " ")
      output$probNoRemission = renderPlot({
        plot(0,type='n',axes=FALSE, main="Welcome to the First Episode Psychosis\nNon-Remission Calculator v0.0.2 (25/11/19)", 
             sub="We used data from the National Evaluation of Development of Early 
          intervention Network study (NEDEN) for model development and 
          internal-external validation. NEDEN is a longitudinal naturalistic 
          study of 1027 patients aged 14-35 with first-episode psychosis
          recruited from 14 early intervention services across the National 
          Health Service (NHS) in England (2005-10).

        For external validation, we used data from the Outlook longitudinal 
        naturalistic study of 399 patients aged 16-35 recruited from an
        additional 11 NHS England early intervention services, throughout 
        April 2006-February 2009."
             , ylab = "",xlab = "",cex.main = 2, col.main = "black")
        text(x = 1, y=0.5, "Developed in collaboration between the Institute of\nMental Health & Wellbeing, University of Glasgow &\nInstitute of Mental Health, University of Birmingham\nby Dr Samuel Leighton and Dr Pavan Mallikarjun.", cex = 1.5)
      })
  })

  observeEvent(input$goButton, {
  output$probNoRemission <- renderPlot({
    #No output before button clicked
    if (input$goButton == 0)
      return()
    
    # Take a dependency on input$goButton
    input$goButton
    
    probNoRemission = isolate({
      
      PCT_Local_Concentration_2007 = getCCGFromPostcode(input$PostCodeToDep, ccg_local_conc_2019)
      
      if(input$Qualification_Level_Ordinal == -1)
      {
        Qualification_Level_Ordinal = NA
      }else
      {
        Qualification_Level_Ordinal = input$Qualification_Level_Ordinal
      }
      
      if(input$Client_Sociability_Withdrawal_Late_Adolescence == -1)
      {
        Client_Sociability_Withdrawal_Late_Adolescence = NA
      }else
      {
        Client_Sociability_Withdrawal_Late_Adolescence = input$Client_Sociability_Withdrawal_Late_Adolescence
      }
      
      if(input$Client_Highest_Functioning_Achieved_Life == -1)
      {
        Client_Highest_Functioning_Achieved_Life = NA
      }else
      {
        Client_Highest_Functioning_Achieved_Life = input$Client_Highest_Functioning_Achieved_Life
      }
      
      if(input$Client_Energy_Level == -1)
      {
        Client_Energy_Level = NA
      }else
      {
        Client_Energy_Level = input$Client_Energy_Level
      }
      
      if(input$BL_PANSS_P3_Hallucinatory_Behaviour == -1)
      {
        BL_PANSS_P3_Hallucinatory_Behaviour = NA
      }else
      {
        BL_PANSS_P3_Hallucinatory_Behaviour = input$BL_PANSS_P3_Hallucinatory_Behaviour
      }
      
      if(input$BL_PANSS_N4_Passive_Social_Withdrawal == -1)
      {
        BL_PANSS_N4_Passive_Social_Withdrawal = NA
      }else
      {
        BL_PANSS_N4_Passive_Social_Withdrawal = input$BL_PANSS_N4_Passive_Social_Withdrawal
      }
      
      if(input$BL_PANSS_G9_Unusual_Thought_Content == -1)
      {
        BL_PANSS_G9_Unusual_Thought_Content = NA
      }else
      {
        BL_PANSS_G9_Unusual_Thought_Content = input$BL_PANSS_G9_Unusual_Thought_Content
      }
      
      if(input$YMRS_Appearance == -1)
      {
        YMRS_Appearance = NA
      }else
      {
        YMRS_Appearance = input$YMRS_Appearance
      }
      
      if(input$IS_Not_Due_Illness == -1)
      {
        IS_Not_Due_Illness = NA
      }else
      {
        IS_Not_Due_Illness = input$IS_Not_Due_Illness
      }
      
      print(paste(Qualification_Level_Ordinal,input$ADJ_DUP,Client_Sociability_Withdrawal_Late_Adolescence,Client_Highest_Functioning_Achieved_Life,
            Client_Energy_Level,BL_PANSS_P3_Hallucinatory_Behaviour,YMRS_Appearance,IS_Not_Due_Illness,input$BL_GAF_Symptoms,PCT_Local_Concentration_2007$DepZScore))
      x = data.frame("Qualification_Level_Ordinal" = as.numeric(Qualification_Level_Ordinal), 
                     "ADJ_DUP" = as.numeric(input$ADJ_DUP), 
                     "Client_Sociability_Withdrawal_Late_Adolescence" = as.numeric(Client_Sociability_Withdrawal_Late_Adolescence),
                     "Client_Highest_Functioning_Achieved_Life" = as.numeric(Client_Highest_Functioning_Achieved_Life),
                     "Client_Energy_Level" = as.numeric(Client_Energy_Level),
                     "BL_PANSS_P3_Hallucinatory_Behaviour" = as.numeric(BL_PANSS_P3_Hallucinatory_Behaviour),
                     "BL_PANSS_N4_Passive_Social_Withdrawal" = as.numeric(BL_PANSS_N4_Passive_Social_Withdrawal),
                     "BL_PANSS_G9_Unusual_Thought_Content" = as.numeric(BL_PANSS_G9_Unusual_Thought_Content),
                     "YMRS_Appearance" = as.numeric(YMRS_Appearance),
                     "IS_Not_Due_Illness" = as.numeric(IS_Not_Due_Illness),
                     "BL_GAF_Symptoms" = as.numeric(input$BL_GAF_Symptoms),
                     "PCT_Local_Concentration_2007" = as.numeric(PCT_Local_Concentration_2007$DepZScore))
      print(str(x))
      tryCatch({
        y = predict.train(mod_eden_all_final_Rem_outlook_glm_new_shrink, x, type = "prob", na.action = na.pass)
        y
      },
      error=function(error_message) {
        print("Cannot impute when all predictors are missing in the new data point")
        NULL
      })
    })
    if(is.null(probNoRemission))
    {
      plot(0,type='n',axes=FALSE, main="Error!\nAll baseline variables are missing!", ylab = "",xlab = "",cex.main = 2, col.main = "red")
    }else
    {
      smileyDiagram(probNoRemission$No[1], estimated = anyNA(x), postCode = PCT_Local_Concentration_2007$completePostCode, CCG = PCT_Local_Concentration_2007$CCGName)
    }
  })
  })
  
}

shinyApp(ui = ui, server = server)
