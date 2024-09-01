const path = require('path');

module.exports = {
  mode: 'development',
  entry:
  {
    main: './webkit/js/main.jsx',

  },
  output: {
    path: path.join(__dirname, '/webkit/static/js/'),
    filename: '[name].bundle.js',
  },
  module: {
    rules: [
      {
        test: /\.(js|jsx)$/,
        exclude: /node_modules/,
        use: {
          loader: "babel-loader"
        }
      },
      {
        test: /\.css$/,
        use: ['style-loader', 'css-loader'],
      },
      {
        test: /\.svg$/,
        use: {
          loader: "file-loader"
        }
      }
    ]
  },
  resolve: {
    extensions: ['.js', '.jsx'],
  },
};
