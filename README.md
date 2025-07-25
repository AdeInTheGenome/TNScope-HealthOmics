 #   S R S   S h o r t R e a d   M a p   &   V a r i a n t   C a l l i n g   ( M u l t i s a m p l e )   -   A W S   H e a l t h O m i c s   W o r k f l o w 
 
 T h i s   A W S   H e a l t h O m i c s   p r i v a t e   w o r k f l o w   ( w r i t t e n   i n   W D L )   p e r f o r m s   r e a d   m a p p i n g   a n d   v a r i a n t   c a l l i n g   o n   p a i r e d   t u m o r - n o r m a l   s h o r t   r e a d   F A S T Q   f i l e s .   I t   s u p p o r t s   m u l t i p l e   s a m p l e s   i n   a   b a t c h   r u n   a n d   u s e s   i n d u s t r y - s t a n d a r d   v a r i a n t   c a l l e r s   ( S e n t i e o n ' s   T N S c o p e   a n d   G o o g l e ' s   D e e p S o m a t i c ) . 
 
 - - - 
 
 # #   O v e r v i e w 
 
 -   * * W o r k f l o w   I D : * *   ` 1 6 5 1 5 6 5 ` 
 -   * * L a t e s t   S u c c e s s f u l   R u n   I D : * *   ` 2 4 2 9 2 1 3 ` 
 -   * * O u t p u t   L o c a t i o n : * *   [ R u n   2 4 2 9 2 1 3   O u t p u t s   i n   S 3 ] ( h t t p s : / / s 3 . c o n s o l e . a w s . a m a z o n . c o m / s 3 / b u c k e t s / s r s - p o c - t e s t / o m i c s - t e s t - a d e / s h o r t r e a d _ p i p e l i n e _ m u l t i s a m p l e _ t e s t / 2 4 2 9 2 1 3 / ? r e g i o n = u s - e a s t - 1 & t a b = o b j e c t s ) 
 -   * * R u n t i m e * *   1 2   h o u r s   a n d   4 0   m i n u t e s   o n   a   t w o   s a m p l e s   ( c h r 4 & c h r 7   i n p u t   f a s t q s ) 
 T h i s   p i p e l i n e : 
 -   A l i g n s   t u m o r   a n d   n o r m a l   F A S T Q   f i l e s   u s i n g   S e n t i e o n   D N A S e q 
 -   C a l l s   v a r i a n t s   u s i n g   T N S c o p e   a n d   D e e p S o m a t i c 
 -   S u p p o r t s   m u l t i s a m p l e   i n p u t   w i t h   s h a r e d   o r   u n i q u e   n o r m a l   s a m p l e s 
 
 - - - 
 
 # #   R e p o s i t o r y   C o n t e n t s 
 
 -   ` S R S - S h o r t R e a d - M a p - V C - M u l t i s a m p l e . w d l `      W o r k f l o w   d e f i n i t i o n 
 -   ` S R S _ S h o r t r e a d _ M a p _ V C _ M u l t i s a m p l e _ p a r a m e t e r s _ d e f i n i t i o n . j s o n `      I n p u t   s c h e m a   a n d   d e s c r i p t i o n s 
 -   ` S R S _ S h o r t r e a d _ M a p _ V C _ M u l t i s a m p l e _ p a r a m e t e r s . j s o n `      S a m p l e   i n p u t   J S O N   f i l e 
 
 - - - 
 
 # #   R e q u i r e d   I n p u t s 
 
 # # #   P a r a m e t e r s 
 
 c o h o r t :   A n   a r r a y   t h a t   d e f i n e s   t u m o r - n o r m a l   F A S T Q   f i l e   g r o u p s   a l o n g   w i t h   s a m p l e   m e t a d a t a .   R e q u i r e d . 
 r e f e r e n c e _ n a m e :   T h e   g e n o m e   r e f e r e n c e   t o   u s e   ( o p t i o n s   i n c l u d e   h g 3 8 ,   t 2 t _ m a t ,   a n d   t 2 t _ p a t ) .   R e q u i r e d . 
 c a n o n i c a l _ u s e r _ i d :   Y o u r   A W S   c a n o n i c a l   u s e r   I D ,   n e e d e d   f o r   o b t a i n i n g   a   S e n t i e o n   l i c e n s e .   R e q u i r e d . 
 s e n t i e o n _ d o c k e r :   U R I   f o r   t h e   S e n t i e o n   D o c k e r   i m a g e   t o   u s e   i n   t h e   w o r k f l o w .   R e q u i r e d . 
 d e e p s o m a t i c _ d o c k e r :   U R I   f o r   t h e   D e e p S o m a t i c   D o c k e r   i m a g e .   R e q u i r e d . 
 n _ t h r e a d s :   N u m b e r   o f   v C P U s   t o   a l l o c a t e   ( d e f a u l t   i s   3 2 ) .   O p t i o n a l . 
 m e m o r y :   A m o u n t   o f   m e m o r y   t o   a l l o c a t e   f o r   t a s k s   ( d e f a u l t   i s   6 4   G i B ) .   O p t i o n a l . 
 
 - - - 
 
 # #   S a m p l e   I n p u t   S t r u c t u r e 
 
 ` ` ` j s o n 
 { 
     " c o h o r t " :   { 
         " s a m p l e s " :   [ 
             { 
                 " t u m o r _ s a m p l e _ n a m e " :   " t e s t _ s a m p l e 1 " , 
                 " n o r m a l _ s a m p l e _ n a m e " :   " c h r 4 _ c h r 7 _ n o r m a l " , 
                 " r 1 _ t u m o r _ f a s t q s " :   [ 
                     " s 3 : / / s r s - p o c - t e s t / f a s t q s / c o m p r e s s e d / c h r 4 _ c h r 7 _ t u m o r _ R 1 . f a s t q . g z " 
                 ] , 
                 " r 2 _ t u m o r _ f a s t q s " :   [ 
                     " s 3 : / / s r s - p o c - t e s t / f a s t q s / c o m p r e s s e d / c h r 4 _ c h r 7 _ t u m o r _ R 2 . f a s t q . g z " 
                 ] , 
                 " r 1 _ n o r m a l _ f a s t q s " :   [ 
                     " s 3 : / / s r s - p o c - t e s t / f a s t q s / c o m p r e s s e d / c h r 4 _ c h r 7 _ n o r m a l _ R 1 . f a s t q . g z " 
                 ] , 
                 " r 2 _ n o r m a l _ f a s t q s " :   [ 
                     " s 3 : / / s r s - p o c - t e s t / f a s t q s / c o m p r e s s e d / c h r 4 _ c h r 7 _ n o r m a l _ R 2 . f a s t q . g z " 
                 ] , 
                 " t u m o r _ r e a d _ g r o u p s " :   [ 
                     " @ R G \ \ t I D : c h r 4 _ c h r 7 _ t u m o r _ t e s t 1 \ \ t S M : c h r 4 _ c h r 7 _ t u m o r _ t e s t 1 \ \ t P L : I L L U M I N A " 
                 ] , 
                 " n o r m a l _ r e a d _ g r o u p s " :   [ 
                     " @ R G \ \ t I D : n o r m a l _ s a m p l e \ \ t S M : c h r 4 _ c h r 7 _ n o r m a l \ \ t P L : I L L U M I N A " 
                 ] 
             } 
         ] 
     } , 
     " r e f e r e n c e _ n a m e " :   " h g 3 8 " , 
     " c a n o n i c a l _ u s e r _ i d " :   " Y O U R _ C A N O N I C A L _ I D " , 
     " s e n t i e o n _ d o c k e r " :   " Y O U R _ S E N T I E O N _ I M A G E _ U R I " , 
     " d e e p s o m a t i c _ d o c k e r " :   " Y O U R _ D E E P S O M A T I C _ I M A G E _ U R I " 
 } 
 ` ` ` 
 
 
 # #   N e e d s 
 
 S e q u e n c i n g   S a m p l e   L o g   t o   b e   u s e d   f o r   g e n e r a t i n g   i n p u t   p a r a m e t e r   f i l e s   s h o u l d   h a v e   c o l u m n s   t h a t   d e s i g n a t e : 
 -   S a m p l e   T y p e   ( T u m o r   o r   N o r m a l ) 
 -   C o r r e s p o n d i n g   N o r m a l   S a m p l e   N a m e   ( N / A   i f   s a m p l e   i n   q u e s t i o n   i s   t h e   n o r m a l ) 
 -   S 3   p a t h s   t o   r 1   f a s t q   a n d   r 2   f a s t q 
