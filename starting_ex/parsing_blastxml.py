#create the dictionary with BLAST results
blast_resdi = {'name':[], 'id':[], 'q_cover':[], 'e_val':[], 'hit_length':[]}

for hit in blast_qresult:
        
    #the hit information
    hitdit = hit[0]
    
    #the hit id cleaned up
        #this might need further cleaning depending on the case
    hitid = str(hitdit.hit_id).replace('ref|', '')
    hitid = hitid.replace('|', '')
    
    #this q_cover metric is: 
        #the length of matching sequence between reference and hit over the length of the reference
    query_cov = (hitdit.query_span/qlen)

    #the length of the matching sequence
    hitlen = len(hitdit.hit) 

    
    #the cutoff here is unnecessary since it can be applied in a later stage
        #if query_cov > 0.5:
      
    
    #add everything to the dict
    
    blast_resdi['name'].append(str(hitdit.hit_description))

    blast_resdi['id'].append(hitid)

    blast_resdi['q_cover'].append(query_cov)
    
    blast_resdi['e_val'].append(hitdit.evalue)

    blast_resdi['hit_length'].append(hitlen)
    