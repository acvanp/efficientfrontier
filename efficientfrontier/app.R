#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#



library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage( shinyjs::useShinyjs(),
              titlePanel("Financial Asset Portfolio Analysis"),
                
              conditionalPanel(condition="input.plot_tabs != 'User guide'",
                                 
              tabsetPanel(id="ui_tab",
              
              tabPanel("Portfolio Optimization", 
                       
                  
                    
                  sidebarPanel(
                  
                  actionButton("update", "Plot Results"),
                  helpText("Plot updates at bottom in mobile browser"),
                  
                  # Copy the line below to make a text input box
                  textInput("tickers", "Tickers list", value = "GOOG, AMZN, TSLA, IVV"),
                  helpText("list up to 8 tickers separated by comma and space like AA, BB, CC"),
                  
                  # Copy the line below to make a date selector 
                  dateInput("start.date", "Info start", value = "2017-01-01"),
                  helpText("information window start date"),
                  
                  # Copy the line below to make a date selector 
                  dateInput("end.date", "Info end", value = "2020-03-15"),
                  helpText("information window end date"),
                  
                  numericInput("inv", "Investment $:", 20000, min = 1500, max = 1000000000),
                  helpText(""),
                  
                  numericInput("risk.percentile", "Risk Percentile (0-100):", 10, min = 1, max = 100),
                  helpText("Percentile out of\ndata point distribution"),
                  
                  # Copy the line below to make a date selector 
                  dateInput("buy.in.date", "Buy in date", value = "2020-03-15"),
                  helpText(""),
                  
                  # Copy the line below to make a date selector 
                  dateInput("cash.out.date", "Cash out date", value = "2020-10-20"),
                  helpText(""),
                  
                  uiOutput("urllinks")
                  
                )
                
                ,
                
                # Show a plot of the generated distribution
                mainPanel(fluidRow(column(8,
                                          h4("Efficient Frontier"),
                                          pre(id = "console"),
                                          plotOutput("efrnt")),
                
                # Show a plot of the generated distribution
                column(8, h4("Recommendation Table:"), tableOutput("tableout")
                
                )
              ))
                  )
              
              ,
              
              tabPanel("Portfolio Tester",
                   
                     
                         sidebarPanel(
                                        
                                        actionButton("update2", "Plot Results"),
                                        helpText("Plot updates at bottom in mobile browser"),
                                        
                                        # Copy the line below to make a text input box
                                        textInput("tickers2", "Tickers list", value = "GOOG, AMZN, TSLA, IVV"),
                                        helpText("list up to 8 tickers separated by comma and space like AA, BB, CC"),
                                        
                                        
                                        # Copy the line below to make a text input box
                                        textInput("shares.count", "Shares bought", value = "3, 3, 3, 3"),
                                        helpText("list the number of shares bought corresponding to each stock ticker listed
                                                 separated by comma and space like 4, 5, 6")
                                        
                                        
                                        #uiOutput("urllinks")
                                        
                                      )
                            
                                      ,
                                      
                                      # Show a plot of the generated distribution
                                      mainPanel(fluidRow(column(8,
                                                                h4("Portfolio Tester"),
                                                                pre(id = "console2"),
                                                                plotOutput("ptester")),
                                      
                                      # Show a plot of the generated distribution
                                      column(8, h4("Portfolio Comparison Table:"), tableOutput("tableout2")
                                      )
                                      )
                                      

                                  )
                   )
              )
  )
)


##############################################

#_________________
### START
#------------------

efficient.frontier = function(tickers, 
                              start.date, 
                              end.date,
                              inv,
                              risk.percentile,
                              buy.in.date,
                              cash.out.date){
  
  
  
  # https://www.codingfinance.com/post/2018-03-27-download-price/
  library(markdown)
  library(quantmod)
  #options("getSymbols.warning4.0"=FALSE)
  #options("getSymbols.yahoo.warning"=FALSE)
  library(gtools)
  library(data.table)
  library(scales)
  library(ggplot2)
  library(stringi)
  
  tickers = unlist(strsplit(tickers, ","))
  tickers = stri_replace_all_fixed(tickers, " ", "")
  tickers = toupper(tickers)
  
  weight.intervals = c(0.01, 0.01, 0.01,
                       0.02, 0.02, 0.1,
                       0.1, 0.1, 0.2,
                       0.2)
  
  n = length(tickers)
  
  x = seq(0,1,weight.intervals[n])
  
  ll = list()
  for(i in 1:(n-1)){
    ll[[i]]  = x  
  }
  
  
  message("Building weights...")
  
  d2 = expand.grid(ll)
  
  #nrow(unique(d2)) == nrow(d2)
  d2$z = 1-rowSums(d2)
  d2 = d2[which(rowSums(d2) <= 1),]
  d2 = d2[which(d2$z  >= 0),]
  
  colnames(d2) = chr(97:122)[1:n]
  
  
  message("Getting stocks...")

  
  getSymbols(tickers,
             from = start.date,
             to = end.date)
  
  prices <- lapply(tickers, function(x) Ad(get(x)))
  

  
  ipo.dates = c()
  
  for(i in 1:n){
    ipo.dates = c(ipo.dates, index(prices[[i]][1]) )
  }
  
  start.date = as.Date(max(ipo.dates))
  
  getSymbols(tickers,
             from = start.date,
             to = end.date)
  
  prices <- lapply(tickers, function(x) Ad(get(x)))
  
  
  
  #reduce the number of rows for assets that tradeon weekends
  ll = sapply(prices, length)
  minrows = prices[[which(min(ll) %in% ll)]]
  
  for(i in 1:n){
    if(i == which(ll == min(ll))){next}
    if( length(prices[[i]]) == length(minrows)){next}
    prices[[i]] = prices[[i]][which(index(prices[[i]]) %in% index(minrows))]
  }
  
  prices = na.omit(prices)
  
  
  #reduce the number of rows for assets that tradeon weekends
  ll = sapply(prices, length)
  minrows = prices[[which(min(ll) %in% ll)]]
  
  for(i in 1:n){
    if(i == which(ll == min(ll))){next}
    if( length(prices[[i]]) == length(minrows)){next}
    prices[[i]] = prices[[i]][which(index(prices[[i]]) %in% index(minrows))]
  }
  
  
  
  
  df = data.frame(unlist(prices[[1]]))
  for(i in 2:n){
    df = cbind(df, prices[[i]])
  }
  
  # index normalization by first value
  for(i in 1:n){df[,i] = df[,i]/(df[,i][1])}
  
  # calculate the necessary values:
  # I) expected returns
  
  returns = c()
  for(i in 1:n){
    returns = c(returns, assign(paste0("er_",colnames(d2)[i]),mean(df[,i])))
  }
  
  
  
  # II) risk (standard deviation) as a risk measure
  
  risk = c()
  for(i in 1:n){
    risk = c(risk , assign(paste0("sd_",colnames(d2)[i]), sd(df[,i])))
    
  }
  
  # III) covariance
  
  covpairs = expand.grid(colnames(d2),colnames(d2))
  
  indexpairs = expand.grid(seq(1,n), seq(1,n))
  
  v = c()
  u = data.frame()
  for(i in 1:nrow(covpairs)){
    if( !paste(unlist(covpairs[i,]), collapse="") %in% v){
      if( !paste(rev(unlist(covpairs[i,])), collapse="") %in% v){
        if(paste(rev(unlist(covpairs[i,])), collapse="") != paste(unlist(covpairs[i,]), collapse="")){
          v = c(v, paste(rev(unlist(covpairs[i,])), collapse="") )
          u = rbind(u, indexpairs[i,])    
        }
      }
    }
  }
  
  covpairs = c()
  for(i in 1:choose(n,2)){
    covpairs = c(covpairs, assign(paste0("cov_",v[i]), cov(df[, u[i,1]], df[,u[i,2]])))
  }
  
  # calculate the expected returns and standard deviations for the 1000 possible portfolios
  
  d2$er_p = 0
  for(i in 1:n){
    d2$er_p = d2$er_p + d2[,i] * returns[i]
  }
  
  
  d2$sd_p = 0
  for(i in 1:n){
    d2$sd_p = d2$sd_p + (d2[,i]^2 * risk[i]^2)
  }
  
  for(i in 1:choose(n,2)){
    d2$sd_p = d2$sd_p + (2 * covpairs[i] * d2[, u[i,1]] * d2[,u[i,2]])
  }
  d2$sd_p = sqrt(d2$sd_p)
  
  # how much volatility is acceptable?
  quantile(d2$sd_p)
  # 110 at the 50% quantile
  
  risk.percentile = risk.percentile / 100 # pick an exotic percentile of risk
  
  diffs = abs(round(d2$sd_p,2) - round(quantile(d2$sd_p, risk.percentile),2))
  
  prop = d2[ which(diffs == min(diffs) ),]
  prop = prop[which(prop$er_p ==  max(prop$er_p)),]
  prop
  
  label.table = data.table(sd = risk, mean = returns)
  label.table = rbind(label.table, data.frame(sd = as.numeric(quantile(d2$sd_p, risk.percentile)), mean = prop$er_p))
  
  message("plotting...")
  
  p = ggplot() +
    geom_point(data = d2, aes(x = sd_p, y = er_p), size = 1, shape = 1) +
    geom_point(data = data.table(sd = risk, mean = returns),
               mapping = aes(x = sd, y = mean), color = "blue", size = 3, shape = 18) +
    geom_point(data = data.table(x = as.numeric(quantile(d2$sd_p, risk.percentile)), y =  prop$er_p),
               mapping = aes(x = x, y = y), color = "red", size = 4, shape = 16) + 
    # Miscellaneous Formatting
    theme_bw() + theme(legend.position = "none") + 
    ggtitle(paste0("Possible Portfolios with ",n, " Risky Assets \nfrom ", start.date, " - ", end.date)) +
    xlab("Volatility") + ylab("Expected Returns") +
    scale_y_continuous(label = percent, limits = c(min(d2$er_p)*0.8, max(d2$er_p) * 1.2)) +
    scale_x_continuous(label = percent, limits = c(min(d2$sd_p)*0.8, max(d2$sd_p) * 1.2)) +
    
    ggrepel::geom_label_repel(data = label.table,
                              mapping = aes(x = sd, y = mean, label = c(tickers, "E.FRNT")), 
                              color = c(rep("blue", n), "red") , size = 3) +
    
    #ggrepel::geom_label_repel(data = data.table(x = as.numeric(quantile(d2$sd_p, risk.percentile)), y =  prop$er_p),
    #                          mapping = aes(x = x, y = y, label = "E.FRNT"), color = "red", size = 2.8) +
    
    geom_segment(
      mapping=aes(
        x = as.numeric(quantile(d2$sd_p, risk.percentile)), 
        y = min(d2$er_p)*0.8, 
        xend = as.numeric(quantile(d2$sd_p, risk.percentile)), 
        yend = prop$er_p, 
        color = "red"), size = 1.1, linetype = 3) +
    geom_segment(
      mapping = aes(
        x = min(d2$sd_p)*0.8, 
        y = prop$er_p, 
        xend = as.numeric(quantile(d2$sd_p, risk.percentile)), 
        yend = prop$er_p, 
        color = "red"), size = 1.1, linetype = 3) 
  
  
  p
  
  buy.in.date = buy.in.date
  cash.out.date = cash.out.date
  
  getSymbols(tickers,
             from = buy.in.date,
             to = cash.out.date)
  
  prices <- lapply(tickers, function(x) Ad(get(x)))
  
  
  
  ipo.dates = c()
  
  for(i in 1:n){
    ipo.dates = c(ipo.dates, index(prices[[i]][1]) )
  }
  
  start.date = as.Date(max(ipo.dates))
  
  getSymbols(tickers,
             from = start.date,
             to = end.date)
  
  prices <- lapply(tickers, function(x) Ad(get(x)))
  
  
  
  #reduce the number of rows for assets that tradeon weekends
  ll = sapply(prices, length)
  minrows = prices[[which(min(ll) %in% ll)]]
  
  for(i in 1:n){
    if(i == which(ll == min(ll))){next}
    if( length(prices[[i]]) == length(minrows)){next}
    prices[[i]] = prices[[i]][which(index(prices[[i]]) %in% index(minrows))]
  }
  
  prices = na.omit(prices)
  
  
  #reduce the number of rows for assets that tradeon weekends
  ll = sapply(prices, length)
  minrows = prices[[which(min(ll) %in% ll)]]
  
  for(i in 1:n){
    if(i == which(ll == min(ll))){next}
    if( length(prices[[i]]) == length(minrows)){next}
    prices[[i]] = prices[[i]][which(index(prices[[i]]) %in% index(minrows))]
  }
  
  
  
  
  df = data.frame(prices[[1]])
  for(i in 2:length(tickers)){
    df = cbind(df, prices[[i]])
  }
  
  
  buyin = df[which(as.Date(rownames(df)) - as.Date(buy.in.date) == min(as.Date(rownames(df)) - as.Date(buy.in.date))),]
  
  cashout = df[which(as.Date(rownames(df)) - as.Date(cash.out.date) == max(as.Date(rownames(df)) - as.Date(cash.out.date))),]
  
  inv = inv # investment $
  
  tickers
  shares.bought = round(inv * prop[,1:n] / buyin)
  
  distribution = round(inv * prop[,1:n] / buyin) * buyin
  
  # percent return
  returns = (cashout - buyin) / buyin
  
  gain = distribution*returns
  
  #print(sum(gain))
  
  pct.returns = round(100*sum(gain)/inv) 
  
  #print(paste(pct.returns, "% return"))
  
  table.out = data.frame( c(tickers, "overall"),
                          c(as.vector(unlist(distribution)), sum(distribution)), 
                          c(as.vector(unlist(shares.bought)), sum(shares.bought)),
                          c(as.vector(round(unlist(buyin),2)), ""), 
                          c(as.vector(round(unlist(cashout),2)), paste0("$", round(sum(gain)), " profit") ), 
                          c(100*round(as.vector(unlist(returns)),2), paste(pct.returns, "% return") ))
  
  colnames(table.out) = c("ticker", "dist.$", "shares", "Buy price $", "Sell price $", "returns %")
  rownames(table.out) = c(tickers, "overall")
  
  
  message("Done!")
  
  return(list(p, table.out))
  
  
}

# 
# # function test area
# ex = efficient.frontier(
#   tickers = c("FB, GOOG, VOO, IVV"), 
#   start.date = "2018-01-01", 
#   end.date = "2020-03-15", 
#   inv = 20000,
#   risk.percentile = 20, 
#   buy.in.date = "2020-03-15", 
#   cash.out.date = "2020-10-20")






###############################################
#####



portfolio.tester = function(tickers, 
                           shares.count){
  
  
  
  # https://www.codingfinance.com/post/2018-03-27-download-price/
  library(markdown)
  library(quantmod)
  #options("getSymbols.warning4.0"=FALSE)
  #options("getSymbols.yahoo.warning"=FALSE)
  library(gtools)
  library(data.table)
  library(scales)
  library(ggplot2)
  library(stringi)
  
  tickers = unlist(strsplit(tickers, ","))
  tickers = stri_replace_all_fixed(tickers, " ", "")
  tickers = toupper(tickers)
  
  weight.intervals = c(0.01, 0.01, 0.01,
                       0.02, 0.02, 0.1,
                       0.1, 0.1, 0.2,
                       0.2)
  
  n = length(tickers)
  
  x = seq(0,1,weight.intervals[n])
  
  ll = list()
  for(i in 1:(n-1)){
    ll[[i]]  = x  
  }
  
  
  message("Building weights...")
  
  d2 = expand.grid(ll)
  
  #nrow(unique(d2)) == nrow(d2)
  d2$z = 1-rowSums(d2)
  d2 = d2[which(rowSums(d2) <= 1),]
  d2 = d2[which(d2$z  >= 0),]
  
  colnames(d2) = chr(97:122)[1:n]
  
  
  message("Getting stocks...")
  
  start.date = "2017-01-01"
  end.date = Sys.Date()
  
  getSymbols(tickers,
             from = start.date,
             to = end.date)
  
  prices <- lapply(tickers, function(x) Ad(get(x)))
  
  ipo.dates = c()
  
  for(i in 1:n){
    ipo.dates = c(ipo.dates, index(prices[[i]][1]) )
  }
  
  start.date = as.Date(max(ipo.dates))
  end.date = Sys.Date()
  
  getSymbols(tickers,
             from = start.date,
             to = end.date)
  
  prices <- lapply(tickers, function(x) Ad(get(x)))
  
  #reduce the number of rows for assets that tradeon weekends
  ll = sapply(prices, length)
  minrows = prices[[which(min(ll) %in% ll)]]
  
  for(i in 1:n){
    if(i == which(ll == min(ll))){next}
    if( length(prices[[i]]) == length(minrows)){next}
    prices[[i]] = prices[[i]][which(index(prices[[i]]) %in% index(minrows))]
  }
  
  prices = na.omit(prices)
  
  
  #reduce the number of rows for assets that tradeon weekends
  ll = sapply(prices, length)
  minrows = prices[[which(min(ll) %in% ll)]]
  
  for(i in 1:n){
    if(i == which(ll == min(ll))){next}
    if( length(prices[[i]]) == length(minrows)){next}
    prices[[i]] = prices[[i]][which(index(prices[[i]]) %in% index(minrows))]
  }
  
  
  
  df = data.frame(unlist(prices[[1]]))
  for(i in 2:n){
    df = cbind(df, prices[[i]])
  }
  
  stock.prices = df[nrow(df),]
  
  
  # index normalization by first value
  for(i in 1:n){df[,i] = df[,i]/(df[,i][1])}
  
  # calculate the necessary values:
  # I) expected returns
  
  returns = c()
  for(i in 1:n){
    returns = c(returns, assign(paste0("er_",colnames(d2)[i]),mean(df[,i])))
  }
  
  
  
  # II) risk (standard deviation) as a risk measure
  
  risk = c()
  for(i in 1:n){
    risk = c(risk , assign(paste0("sd_",colnames(d2)[i]), sd(df[,i])))
    
  }
  
  # III) covariance
  
  covpairs = expand.grid(colnames(d2),colnames(d2))
  
  indexpairs = expand.grid(seq(1,n), seq(1,n))
  
  v = c()
  u = data.frame()
  for(i in 1:nrow(covpairs)){
    if( !paste(unlist(covpairs[i,]), collapse="") %in% v){
      if( !paste(rev(unlist(covpairs[i,])), collapse="") %in% v){
        if(paste(rev(unlist(covpairs[i,])), collapse="") != paste(unlist(covpairs[i,]), collapse="")){
          v = c(v, paste(rev(unlist(covpairs[i,])), collapse="") )
          u = rbind(u, indexpairs[i,])    
        }
      }
    }
  }
  
  covpairs = c()
  for(i in 1:choose(n,2)){
    covpairs = c(covpairs, assign(paste0("cov_",v[i]), cov(df[, u[i,1]], df[,u[i,2]])))
  }
  
  # calculate the expected returns and standard deviations for the 1000 possible portfolios
  
  d2$er_p = 0
  for(i in 1:n){
    d2$er_p = d2$er_p + d2[,i] * returns[i]
  }
  
  
  d2$sd_p = 0
  for(i in 1:n){
    d2$sd_p = d2$sd_p + (d2[,i]^2 * risk[i]^2)
  }
  
  for(i in 1:choose(n,2)){
    d2$sd_p = d2$sd_p + (2 * covpairs[i] * d2[, u[i,1]] * d2[,u[i,2]])
  }
  d2$sd_p = sqrt(d2$sd_p)
  
  shares.count = unlist(strsplit(shares.count, ","))
  shares.count = stri_replace_all_fixed(shares.count, " ", "")
  shares.count = as.numeric(shares.count)
  
  distribution.dollars = as.vector(unlist(stock.prices)) *  shares.count
  
  total.distribution = sum(distribution.dollars)
  
  portfolio.weights = distribution.dollars / total.distribution
  
  portfolio.risk = sum(portfolio.weights * risk)
  
  portfolio.er = sum(portfolio.weights * returns)
  
  #risk.percentile = risk.percentile / 100 # pick an exotic percentile of risk
  
  diffs = abs(round(d2$sd_p,2) - portfolio.risk)
  
  prop = d2[ which(diffs == min(diffs) ),]
  prop = prop[which(prop$er_p ==  max(prop$er_p)),]
  prop
  
  label.table = data.table(sd = risk, mean = returns)
  label.table = rbind(label.table, 
                      data.frame(sd = portfolio.risk, mean = prop$er_p),
                      data.frame(sd = portfolio.risk, mean = portfolio.er))
  
  message("plotting...")
  
  p = ggplot() +
    geom_point(data = d2, aes(x = sd_p, y = er_p), size = 1, shape = 1) +
    geom_point(data = data.table(sd = risk, mean = returns),
               mapping = aes(x = sd, y = mean), color = "blue", size = 3, shape = 18) +
    geom_point(data = data.table(x = portfolio.risk, y =  prop$er_p),
               mapping = aes(x = x, y = y), color = "red", size = 2, shape = 0) +
    geom_point(data = data.table(x = portfolio.risk, y =  portfolio.er),
               mapping = aes(x = x, y = y), color = "yellow", size = 2, shape = 2) +
    # Miscellaneous Formatting
    theme_bw() + theme(legend.position = "none") + 
    ggtitle(paste0("Possible Portfolios with ",n, " Risky Assets \nfrom ", start.date, " - ", end.date)) +
    xlab("Volatility") + ylab("Expected Returns") +
    scale_y_continuous(label = percent, limits = c(min(d2$er_p)*0.8, max(d2$er_p) * 1.2)) +
    scale_x_continuous(label = percent, limits = c(min(d2$sd_p)*0.8, max(d2$sd_p) * 1.2)) +
    
    ggrepel::geom_label_repel(data = label.table,
                              mapping = aes(x = sd, y = mean, label = c(tickers, "E.FRNT", "You")), 
                              color = c(rep("blue", n), "red", "black") , size = 3) +
    
    #ggrepel::geom_label_repel(data = data.table(x = as.numeric(quantile(d2$sd_p, risk.percentile)), y =  prop$er_p),
    #                          mapping = aes(x = x, y = y, label = "E.FRNT"), color = "red", size = 2.8) +
    
    geom_segment(
      mapping=aes(
        x = portfolio.risk, 
        y = min(d2$er_p)*0.8, 
        xend = portfolio.risk, 
        yend = prop$er_p, 
        color = "red"), size = 1.1, linetype = 3) +
    geom_segment(
      mapping = aes(
        x = min(d2$sd_p)*0.8, 
        y = prop$er_p, 
        xend = portfolio.risk, 
        yend = prop$er_p, 
        color = "red"), size = 1.1, linetype = 3) 
  
  p
  
  risk.percentiles = c()
  
  for(i in 1:100){
    risk.percentiles = c(risk.percentiles, quantile(d2$sd_p, i/100))
  }
  
  diffs = abs(risk.percentiles - portfolio.risk)
  
  risk.percentile = risk.percentiles[which(diffs == min(diffs))]
  
  risk.percentile = stri_replace_all_fixed(labels(risk.percentile), "%", "")
  
  
  table.out = data.frame(
    c(tickers,
      "Risk percentile", 
      "%Volatility", 
      "Your expected returns %", 
      "Optimal expected returns %"),
    c(shares.count, 
      as.numeric(risk.percentile), 
      100*(round(portfolio.risk, 2)), 
      100*round(portfolio.er,2), 
      100*round(prop$er_p, 2) ), 
    c( round(as.vector(unlist(stock.prices))), rep("", 4)) , 
    c( round(distribution.dollars), rep("", 4)), 
    c(100*round(portfolio.weights,3), rep("", 4)) ,
    c(100*round(as.vector(unlist(prop[1:n])),2), rep("", 4)) 
  )
  
  colnames(table.out) = c("tickers", "shares", "buy price $", "money in $", "your portfolio\nweights", "optimal portfolio\nweights" )
  
  message("Done!")
  
  return(list(p, table.out))
  
  
}



#####
##############################################

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  v <- reactiveValues(data=NULL)
  
  observeEvent(
    
    input$update,
    
    {
      
      
      v$data1 = withCallingHandlers(
        
        
        efficient.frontier(input$tickers, 
                           input$start.date, 
                           input$end.date, 
                           input$inv,
                           input$risk.percentile,
                           input$buy.in.date,
                           input$cash.out.date
        ),
        
        # can use "warning" instead/on top of "message" to catch warnings too
        message = function(m) {
          shinyjs::html("console", m$message, FALSE)
        }
        
      )
      
    }
  )
  
  
  output$efrnt = renderPlot({
    
    v$data1[[1]]
    
  })
  
  output$tableout <- renderTable(v$data1[[2]])
  
  url <- a("GitHub and Info", href = "https://github.com/acvanp/efficientfrontier")
  output$urllinks <- renderUI({
    tagList(url)
  })
  
  
  
  
  u <- reactiveValues(data=NULL)
  
  observeEvent(
    
    input$update2,
    
    {
      
      
      u$data2 = withCallingHandlers(
        
        
        portfolio.tester(input$tickers2, 
                           input$shares.count
        ),
        
        # can use "warning" instead/on top of "message" to catch warnings too
        message = function(m) {
          shinyjs::html("console2", m$message, FALSE)
        }
        
      )
      
    }
  )
  
  
  output$ptester = renderPlot({
    
    u$data2[[1]]
    
  })
  
  output$tableout2 <- renderTable(u$data2[[2]])
  
  
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
