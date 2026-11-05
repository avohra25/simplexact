\ 0/1 knapsack (a 10-item instance, capacity 165).
maximize
  value: 23 a + 31 b + 29 c + 44 d + 53 e + 38 f + 63 g + 85 h + 89 i + 82 j
subject to
  weight: 92 a + 57 b + 49 c + 68 d + 60 e + 43 f + 67 g + 84 h + 87 i + 72 j
          <= 165
binary
  a b c d e f g h i j
end
