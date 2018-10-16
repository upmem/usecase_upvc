#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "upvc.h"

/* 
	Most of this code is inspired by geeksforgeeks.org
*/

int max(int a, int b) 
{
  return (a > b) ? a : b;
}

VARTREE* initialize(int pos, VARIANT* vars) 
{
  VARTREE* root = malloc(sizeof(VARTREE));
  root->pos = pos;
  root->vars = vars;
  root->right = NULL;
  root->left = NULL;
  root->height = 1;
  return root;
}

int height(VARTREE* root) 
{
  if (root == NULL)
    return 0;
  else
    return root->height;
}

int getBalance(VARTREE* root) 
{
  if (root == NULL)
    return 0;
  else
    return height(root->left) - height(root->right);
}

// A utility function to right rotate subtree rooted with y
// See the diagram given above.
VARTREE *rightRotate(VARTREE *y)
{
  VARTREE *x = y->left;
  VARTREE *T2 = x->right;
  
  // Perform rotation
  x->right = y;
  y->left = T2;
  
  // Update heights
  y->height = max(height(y->left), height(y->right))+1;
  x->height = max(height(x->left), height(x->right))+1;
  
  // Return new root
  return x;
}
 
// A utility function to left rotate subtree rooted with x
// See the diagram given above.
VARTREE *leftRotate(VARTREE *x)
{
  VARTREE *y = x->right;
  VARTREE *T2 = y->left;
  
  // Perform rotation
  y->left = x;
  x->right = T2;
  
  //  Update heights
  x->height = max(height(x->left), height(x->right))+1;
  y->height = max(height(y->left), height(y->right))+1;
  
  // Return new root
  return y;
}

VARTREE* insert(VARTREE *root, VARIANT* var) 
{
  // if the tree does not exist, initialize it with the given variant
  if (root == NULL) {
    return initialize(var->pos, var);
  }

  // if the tree is rooted at a variant with the same position, increase the variant depth
  if (root->pos == var->pos)
    {
      root->vars->depth += 1;
      free(var);
      return root;
    }

  // if the variant is at a higher position than the root, then insert it in the right subtree
  if (var->pos > root->pos)
    {	
      root->right = insert(root->right, var);
    }
  // if the variant is at a lower position than the root, then insert it in the left subtree
  else 
    {
      root->left = insert(root->left, var);
    }

  root->height = 1 + max(height(root->left), height(root->right));

  int balance = getBalance(root);

  // If balance magnitude is greater than 1, then the tree needs to be rebalanced

  // Left Left Case
  if (balance > 1 && var->pos < root->left->pos)
    return rightRotate(root);
 
  // Right Right Case
  if (balance < -1 && var->pos > root->right->pos)
    return leftRotate(root);
 
  // Left Right Case
  if (balance > 1 && var->pos > root->left->pos)
    {
      root->left =  leftRotate(root->left);
      return rightRotate(root);
    }
 
  // Right Left Case
  if (balance < -1 && var->pos < root->right->pos)
    {
      root->right = rightRotate(root->right);
      return leftRotate(root);
    }

  // Tree is still balanced
  return root;
}

void insertVariants(LIST_VARIANTS *LV, VARIANT* var)
{
  LV->vt = insert(LV->vt, var);
}

void visuTree(VARTREE *root)
{
  VARIANT *v;

  if (root == NULL)
    return ;

  visuTree(root->left);

  v = root->vars;
  if ((v->depth > 7) && (strlen(v->ref)>1))
  printf("%s\t%d\t.\t%s\t%s\t.\t.\tDEPTH=%d;\n", v->chr, v->pos+1-v->offset, v->ref, v->alt, v->depth);

  visuTree(root->right);

}

void visuVariants(LIST_VARIANTS *LV)
{
  visuTree(LV->vt);
}
