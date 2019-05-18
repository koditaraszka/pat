  #how to compute log mvalues                                                   
  @staticmethod                                                                 
  def sumlog(x, y):                                                             
    print('x: ' + str(x))                                                       
    print('y: ' + str(y))                                                       
    combine = np.stack((x,y), axis=1)                                           
    print('combine: ' + str(combine))                                           
    high = np.amax(combine, axis=1)                                             
    print('high: ' + str(high))                                                 
    low = np.amin(combine, axis=1)                                              
    print('low: ' + str(low))                                                   
    diff = np.subtract(high,low)                                                
    print('diff: ' + str(diff))                                                 
    total = low + np.log1p(np.exp(diff))                                        
    print('total: ' + str(total))                                               
    return total


  #set_mvalues returns the posterior marginal for each trait                    
  def set_mvalues(self,sims):                                                   
    print('sims shape: ' + str(sims.shape))                                     
    all_configs = []                                                            
    for i in range(0,len(self.alt_configs)):                                    
      #print('i: ' + str(i))                                                    
      #print(self.alt_configs[i])                                               
      value = multivariate_normal.logpdf(sims, self.mean, self.alt_configs[i])  
      #value = multivariate_normal.pdf(sims, self.mean, self.alt_configs[i])    
      #print(value[:10])                                                        
      if len(all_configs) is 0:                                                 
        all_configs = value                                                     
      else:                                                                     
        all_configs = Pat.sumlog(all_configs, value)                            
        #all_configs = np.add(all_configs,value)                                
      for j in range(0,len(self.groups)):                                       
        if i in self.groups[j]:                                                 
          if len(self.mvalues[j]) is 0:                                         
            self.mvalues[j] = value                                             
          else:                                                                 
            self.mvalues[j] = Pat.sumlog(self.mvalues[j], value)                
            #self.mvalues[j] = np.add(self.mvalues[j],value)                    
    #print(all_configs[:10])                                                    
    self.mvalues = np.exp(np.subtract(self.mvalues,all_configs))                
    #for k in range(0,len(self.mvalues)):                                       
        #print(self.mvalues[k][:10])                                            
        #self.mvalues[k] = np.exp(self.mvalues[k]-all_configs))                 
    return self.mvalues
