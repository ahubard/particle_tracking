function [s,w]=smbcp(cn,sd,sf,dd,df,u,p)

  [s,w]=unix(sprintf('smbclient //%s/%s "%s" -U %s -c "cd %s;put %s %s"',cn,sd,p,u,dd,sf,df));
