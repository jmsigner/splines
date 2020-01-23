#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(fda)


# Define UI for application that draws a histogram
ui <- fluidPage(
    # Application title
    titlePanel("Splines"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            helpText("Code developed by Roland Langrock (http://www.rolandlangrock.com/spline_smooth.txt)"),
            sliderInput("nb",
                        "Number of basis functions:",
                        min = 1,
                        max = 50,
                        value = 25),
            sliderInput("lambda",
                        "log10(lambda) (penalty):",
                        min = -10,
                        max = 10, step = 0.1,
                        value = 0),
            sliderInput("order",
                        "Order:",
                        min = 1,
                        max = 5,
                        value = 3), 
            helpText("Fork code on github: https://github.com/jmsigner/splines")
        ),

        # Show a plot of the generated distribution
        mainPanel(
            textOutput("test"),
           plotOutput("base_plot"), 
           plotOutput("weighted_base"), 
           plotOutput("reg")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, server) {
    
    cbPalette <- c("#999999","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")
    # simulate data
    set.seed(123)
    x <- runif(100,1,10)
    y <- rnorm(100,sin(1.7*x)+0.5*cos(1.3*x)+0.4*x-0.3*(x-5)^2,0.5)
    # wiggliness penalty
    penalty <- function(as){
        vec <- diff(diff(as))
        return(sum(vec^2))
    }
    # minus penalised sum of squares
    mpls<-function(parvect, lambda, nb, Bx, y){
        gamma <- parvect[1:nb]
        return(sum((y - Bx %*% gamma)^2) + lambda * penalty(gamma))
    }
    
    lambda <- reactive({10^input$lambda})
    
    spline <- reactive({
         degree = input$order - 1 # B-spline degree
         nrknots = input$nb - (input$order - 2) # no. of knots
         knots = seq(-1, 12, length = nrknots + 2 * degree)
         # set up B-spline design matrix
         Bx <- spline.des(knots, x, input$order, outer.ok=TRUE)$design
         mod <- nlm(mpls,rep(1, input$nb), lambda = lambda(), nb = input$nb, Bx = Bx, y = y)
         z <- seq(-1, 11.5, length = 1000)
         B <- spline.des(knots, z, input$order, outer.ok=TRUE)$design
         list(z = z, B = B, nb = input$nb, mod = mod)
    })
    
    
    output$base_plot <- renderPlot({
        par(mfrow=c(2,2),mar=c(2,2,2,2))
        plot(x,y,main="basis functions",pch=19,bty="n")
        for (k in 1:spline()$nb){
            lines(spline()$z, spline()$B[,k],col=cbPalette[k%%7+1],lwd=2)
        }
        
        plot(x, y, main="weighted basis functions", pch=19, bty="n")
        for (k in 1:spline()$nb){
            lines(spline()$z, spline()$mod$estim[k] * spline()$B[,k],
                  col = cbPalette[k%%7+1], lwd = 2)
        }
        
        plot(x,y,main="fitted regression function",pch=19,bty="n")
        lines(spline()$z, spline()$B %*% spline()$mod$estim[1:spline()$nb],
              type="l",main="",lwd=2,col=cbPalette[3])
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
