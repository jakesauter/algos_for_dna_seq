#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Implement the naive_with_rc function
from naive_with_rc import naive_with_rc


# ### Example 1

# In[2]:



p = 'CCC'
ten_as = 'AAAAAAAAAA'
t = ten_as + 'CCC' + ten_as + 'GGG' + ten_as
print(naive(p, t))
print(naive_with_rc(p, t))


# ### Example 2

# In[3]:


p = 'CGCG'
t = ten_as + 'CGCG' + ten_as + 'CGCG' + ten_as
print(naive(p, t))
print(naive_with_rc(p, t))



# ### Example 3

# In[4]:


# Phi-X genome
phix_genome = read_genome('data/phix.fa')


# In[7]:


occurrences = naive_with_rc('ATTA', phix_genome)


# In[8]:


print('offset of leftmost occurrence: %d' % min(occurrences))


# In[9]:


print('# occurrences: %d' % len(occurrences))

