notify_me <- function(condition = TRUE, message = NULL){
  hook_to_me <- "https://hooks.slack.com/services/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"
  
  
  if(is.null(message)){
    if(condition){
      msg <- list(text="Your process has finished!")
      httr::POST(hook_to_me, encode = "json", body = msg)
    }else{
      msg <- list(text="Your process failed :(")
      httr::POST(hook_to_me, encode = "json", body = msg)
    }
  }else{
    if(condition){
      msg <- list(text=message)
      httr::POST(hook_to_me, encode = "json", body = msg)
    }else{
      msg <- list(text="Your process failed :(")
      httr::POST(hook_to_me, encode = "json", body = msg)
    }
  }
  
}